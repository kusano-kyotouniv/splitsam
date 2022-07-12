$| = 1;

# hisat2 や minimap2 で作った .sam を一定数のリファレンスごとに分ける。
# usage: $ perl splitsam.pl mapping.sam 7 150
# 引数１は マッピング結果のファイル。必須。
# 引数２は 何個に分けるか。必須。結果のファイル名は splitsam_0.sam のような感じになる。
# 引数３は 一度に走らせるプロセスの数。メモリが足りない場合は少なくする。省略すると100になる。
# 元のマッピング結果と同サイズ程度の空きHDD容量が必要。マップ率が低い場合はその分HDD要求量が減る。
# 
# samtools が機能しないほどの巨大なデータに使うためのもの。
# ロングリードデータにショートリードデータをマッピングしたような、
# リファレンスもリードも多い巨大なデータなどに。
# 染色体までアセンブルできているような、リファレンスが数本しかない場合は
# 素直にsamtoolsを使ったほうがいい。
# 
# データが大きいので、メモリとCPUをマックス頑張らせてちょっと早くする。
# しかしメモリの搭載量を認識していないので、人力で引数３を調節して動く水準に整える。
# なお、引数３は計算結果には影響しない。少ないと遅い。多いと早い。多すぎるとエラーで止まる。
#
# ７個に分ける場合、７x100 = 700 個のプロセスを１サイクルで走らせる。
# 引数３が150の場合、この700のうち150個が走り出したら次のプロセスは走り出さずに待つ。
# メモリ使用率に余裕がありそうな場合は同時プロセス数はもうすこし増やしたほうが早くなる。
# 
# 動作の例）リファレンス1,387,127個を含む mapping.sam を 引数２=7, （引数３=150） で処理した場合
# ひとつあたり 1,387,127 / 7 / 100 = リファレンス1982 = $ref_unit 個分のファイルに一旦分割する。
# 最後に100個ごとにマージして７個にする。この100は内部で適当に定めた変数。
# 
# 内部で定数にしてるのはこれ。$data_unit は引数にしといた方が便利だろうか。
# $splitfiles_per_merge = 100	一時ファイルの数はこれと引数２の分割数を掛けた値。多めに分ける方が早いだろうか？
# $data_unit = 5000000	一度に処理するデータ量。メモリが足りなくてCPUがあまるときはこれを減らす。
#
# run_polish.sh は hisat2 -> splitsam -> hypo -> (merge) を行うもの。
# 十分な容量のあるハードディスクを計算用紙として指定して使うため、パスを入力する。
# リファレンスデータをsplitsamの出力と同じ形で分割しておく必要がある。
# split_ontreads.pl に、splitsam.plのレポートに現れる $ref_unit の値（1982とか）を入れて実行する。

my $mapfile = $ARGV[0];		# 引数１はマップファイルの名前
if( -f $mapfile){
	print 'mapfile: '.$mapfile."\n";;
	my $filesize = -s $mapfile;
	$filesize = int($filesize / 1000000000);
	print 'mapfile size: '.$filesize." Gb\n";
}else{
	print 'mapfile not found'."\n";
	exit();
}

my $argv_outfile = $ARGV[1];		# 引数２は最終出力ファイルの数
my $outfile_merge_prefix = 'splitsam_';
if($argv_outfile > 1){
	print 'output files:'."\n";
	for(my $i=0;$i<$argv_outfile;$i++){
		print $outfile_merge_prefix.$i.'.sam'."\n";
	}
}else{
	print 'please input number of output files (ex. 6) in second argument'."\n";
	exit();
}
my $merged_outfiles = $argv_outfile;

my $argv_process = $ARGV[2];		# 引数３は同時進行するプロセスの数
my $simultaneous_processes = 100;	# 同時進行するプロセス数のリミット。
if($argv_process > 1){
	$simultaneous_processes = $argv_process;
}else{
	print 'please input threads in multiprocess (ex. 100) in third argument'."\n";
	exit();
}

open IN, $mapfile;			# リファレンスのリード数を数える
my $mapfile_refs = 0;
RCN:while(my $line = <IN>){
	my @s = split(/\t/,$line);
	if($s[0] eq '@SQ'){$mapfile_refs++;}
	if($s[0] eq '@PG'){
		print 'reference seqs: '.$mapfile_refs."\n";
		last RCN;
	}
	if($mapfile_refs % 100000 == 0){print 'counting refs: '.$mapfile_refs."\n";}
}
close IN;
my $splitfiles_per_merge = 100;	# splitいくつ単位でmergeするか。
my $ref_unit = int($mapfile_refs / $merged_outfiles / $splitfiles_per_merge)+1;

open LOG, '>>splitsam_log.txt';
	print LOG 'ARGV1: mapfile: '."\t".$mapfile."\n";
	print LOG 'ARGV2: output files:  '."\t".$merged_outfiles."\n";
	print LOG 'ARGV3: processes:     '."\t".$simultaneous_processes."\n";
	print LOG 'reference seqs:       '."\t".$mapfile_refs."\n";
	print LOG 'refs per internalfile:'."\t".$ref_unit."\n";
	print LOG 'files per finalout:   '."\t".$splitfiles_per_merge."\n";
close LOG;


my $outfile_prefix = '_split';	# 出力ファイルの名前。これに番号がつく。
my $refname_length = 36 + 2;	# リファレンス名の名前の文字数。porechop で２文字増えるので、比較に使うのは38文字+区切り文字1 = 39 文字になる。
my $refstr_length = $refname_length+1;	# 区切り文字の空白1文字分足す。


# リファレンス部分。IDを並べて長い文字列にする。いくつかに分けて作る。
open IN, $mapfile;
my $hd = <IN>;	# 最初の行はヘッダー。最初に１行だけなのでこれでOK
my $pg = '';	# コレはリファレンス部分の最後に出てくるのでそのときロードする。
my @ref_searchid;	# $ref_searchid[$aliquots]	# IDを並べた文字列。
my @refline;		# $refline[$aliquots][$refs]		# リファレンス行を保存しておく。出力用。
my $refs=0;		# リファレンスの数
my @aliquot_size;	# $aliquot_size[$aliquots]	# 各単位に含まれるリファレンスの数。最後だけ $unit と違う数になる
my $aliquots=0;	# 分割数。
my $ref_loop_flag = 1;
while($ref_loop_flag == 1){
	$aliquot_size[$aliquots]=0;
	my @refid;
	SQ:for(my $i=0;$i<$ref_unit;$i++){
		my $line =<IN>;
		my @s = split(/\t/,$line);
		if($s[0] eq '@PG'){$pg = $line;$ref_loop_flag=0;last SQ;}
		$refid[$i] = substr($s[1],3,$refname_length);
		if(length($refid[$i]) == $refname_length-2){$refid[$i] = $refid[$i].'__';}
		$refline[$aliquots][$aliquot_size[$aliquots]] = $line;
		$aliquot_size[$aliquots]++;
		$refs++;
		# if($refs % 2000 ==0){print '.';}
	}
	$ref_searchid[$aliquots] = join('|',@refid);
	@refid = ();
	print 'aliquot #'.$aliquots.' sized '.$aliquot_size[$aliquots]."\n";
	$aliquots++;
}

# 一時的出力先ファイルを生成して、各リファレンス部分を出力しておく仕様だったのを変更。
# 一時的出力先ファイルにはヘッダー部分が無い。
for(my $a=0;$a<$aliquots;$a++){
	open OUT, '>'.$outfile_prefix.$a.'.sam';
#	print OUT $hd;
#	for(my $r=0;$r<$aliquot_size[$a];$r++){print OUT $refline[$a][$r];}
#	print OUT $pg;
	close OUT;
}
print 'created '.$aliquots.' .sam files with prefix '.$outfile_prefix."\n";

# merge後のファイルを生成してヘッダー部分だけ作っておく
for(my $m=0;$m<$merged_outfiles;$m++){
	open OUT, '>'.$outfile_merge_prefix.$m.'.sam';
	print OUT $hd;
	my $mm = $splitfiles_per_merge * $m;
	for(my $a=0;$a<$splitfiles_per_merge;$a++){
		for(my $r=0;$r<$aliquot_size[$a+$mm];$r++){
			print OUT $refline[$a+$mm][$r];
		}
	}
	print OUT $pg;
	close OUT;
}

undef @refline;
undef @aliquot_size;

# 続きはデータ部分。適当な単位でデータをHDDからロードする。

my $data_unit =  5000000;
my $data_loop_flag=1;
my @data;		# data[$datas]		# $datas は unmapped をスキップするので $data_unit より小さくなる。
my @data_id;		# data_id[$datas]	# 判定に使う部分だけ抜いたもの
my $data_id_amphora;	# data_id_amphora[$datas/100]	# or演算子で並べてみる
my $loops=0;
my $stdout_freq = int($data_unit/100);
while($data_loop_flag == 1){	# 一定の単位数のデータをロード>照合>セーブする、のループ。データの終わりまで行ったら終わり。
	print "\n".'starting loop #'.$loops.' calculating '.$data_unit.' lines'."\n";
	# データをロードする。
	my $datas=0;
	my $lines=0;
	DAT:for($d=0;$d<$data_unit;$d++){
		if(my $line =<IN>){
			my @s = split(/\t/,$line);
			if($s[2] ne '*'){
				$data[$datas] = $line;
				$data_id[$datas] = substr($s[2],0,$refname_length);
				if(length($data_id[$datas]) == $refname_length-2){
					$data_id[$datas] = $data_id[$datas].'__';
				}
				$datas++;
			}
		}else{$data_loop_flag = 0; last DAT;}
		if($lines % $stdout_freq == 0){print '.';}
		$lines++;
	}
	print $datas.' mapped reads loaded from '.$lines."\n";
	
	# 検索に使うデータのヒット先refのIDをならべる。こっちが検索対象になる。区切りは空白。
	$data_id_amphora = join(' ',@data_id);
	@data_id = ();
	
	# 並列計算できるか？
	my @pid;
	my $pids=0;
	my $process_count = 0;
	for(my $t=0;$t<$aliquots;$t++){
		$pid[$t] = fork;
		if($pid[$t] == 0){
			$ret = &matchout($t);
			print 'loop #'.$loops.' fork #'.$t.' finish: '.$ret.' hits /'.$datas."\n";
			exit();
		}else{
			$pids++;
			$process_count++;
		}
		while($process_count > $simultaneous_processes){
			my $str = wait;
			$process_count--;
			# print 'Memory space allocation #'.$str.'　released: '.$process_count.' / '.$pids.' processes waiting'."\n";
		}
	}
	
	while($process_count > 0){	# 全部のforkを走らせ終わったら、全部が終わるまで待つ。
		my $str = wait;
		$process_count--;
		print 'Memory space allocation #'.$str.'　released: '.$process_count.' / '.$pids.' processes waiting'."\n";
	}
	
	@data = ();		# 変数自体はまた使うので、メモリの解放に undef しない。
	$data_id_amphora = '';
	
	$loops++;
}
close IN;
print 'split into internal files normally finish'."\n";


# 一旦分割したファイルを指定の個数にマージする。ヘッダーはすでに出力済みなのでリードデータ部分だけ。
for(my $m=0;$m<$merged_outfiles;$m++){

	my $start_num = $splitfiles_per_merge * $m;
	my $end_num = $start_num + $splitfiles_per_merge;
	
	open OUT, '>>'.$outfile_merge_prefix.$m.'.sam';
	# 本体のデータを出力
	for(my $f=$start_num;$f<$end_num;$f++){
		print $outfile_prefix.$f.'.sam';
		open IN, $outfile_prefix.$f.'.sam';
		my $body = do{local $/; <IN>};
		print OUT $body;
		close IN;
		print ' done'."\n";
		unlink($outfile_prefix.$f.'.sam');
	}
	close OUT;
}

exit();


# マルチスレッドで呼ばれる関数。
sub matchout{
	my $id = $_[0];
	my $hits=0;
	my @hit;
	while($data_id_amphora =~ m/$ref_searchid[$id]/g){
		$hit[$hits] = (pos($data_id_amphora) +1) / $refstr_length -1;
		$hits++;
	}
	open OUT, '>>'.$outfile_prefix.$id.'.sam';
	for(my $h=0;$h<$hits;$h++){
		print OUT $data[$hit[$h]];
	}
	close OUT;
	undef @hit;
	return($hits);
}



