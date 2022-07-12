# Nanopore リードデータを分ける。ポリッシングのサイクルの初期データ。
# 分ける元はendtrimとporechopで処理したもの。マッピングツールにも与えられるfasta１ファイルになっている
# splitsam と使うときは定数をそっちにあわせる。

open IN, '/media/k/HDD5/220617yewgenome_polished.fasta';
my $outfiles = 7;			# 出力するファイルの数。		600個に分けて60個ずつまとめて10個にしてうｒ
my $splitmaps = 700;			# splitsam9 でマップを分けた数	600個に分けてる
#my $splitmaps_per_merge = $splitmaps/$outfiles;	# 出力ファイル１個あたりのONTリード数=splitmap数
my $splitmaps_per_merge = 100;		# 手動で数字いれたほうが正確そうなので。100 splitmap ずつ。	
my $refs_per_splitmap = 1982;		# splitmap　ひとつあたりに含まれるrefの数

my $ontreads_per_file = $refs_per_splitmap * $splitmaps_per_merge;	# 出力ファイルひとつあたりのリード数。
for(my $r=0;$r<$outfiles;$r++){
	open OUT, '>splithypo_'.$r.'.fa';
	my $count=0;
	LOOP:for(my $i=0;$i<$ontreads_per_file;$i++){
		my $name=<IN>;
		my $seq=<IN>;
		if($name =~ /\>/){
			print OUT $name;
			print OUT $seq;
			$count++;
		}else{
			last LOOP;
		}
	}
	print '>splithypo_'.$r.'.fa'."\t".$count.' ONTreads'."\n";
	close OUT;
}
close IN;

