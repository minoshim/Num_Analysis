# これは何ですか？
ここは、[計算物理春の学校2026](https://compphysschool.github.io/2026/) 個別講義D-4 「無衝突系の運動論シミュレーション入門」で用いるプログラム置き場です。<br>
- FDS: 有限差分法
- FVS: 有限体積法
- Vlasov: 運動論シミュレーション
- common: 共通C++プログラムとPython可視化補助

C++とPythonプログラムが準備されていますので、お好きな方をご利用ください。<br>
そのまま使うもよし、機能を追加するもよし、自作の参考にするもよしです。<br>
講義ではPythonを用います。<br>

# どうやって使いますか?
ローカル環境で使う場合は、ご自身のPCにファイルをダウンロードするか、
```
git clone https://github.com/minoshim/Num_Analysis.git
```
してください。<br>
クラウド環境で使う場合は、Google Colabを使いますので、Googleアカウントをご準備ください。
## Python
ディレクトリ`Py`の先にある`main.py`がバッチファイル、`hoge.ipynb`がノートブックです（一部、`main.py`のみの補助プログラムもあります）。<br>
バッチファイルをローカル環境で使う場合は、numpy, scipy, matplotlibがインストール済みのPythonを準備します。<br>
使いたい課題のディレクトリに移動してから、例えば
```
cd FDS/Py/Adv
python main.py
```
のように実行します。Pythonの対話環境から実行する場合は、
```
exec(open("main.py").read())
```
とすると計算が行われます。

ノートブックをクラウド環境で使う場合は、例えば
```
https://github.com/minoshim/Num_Analysis/blob/main/FDS/Py/Adv/fds_adv.ipynb
```
のURLを
```
https://colab.research.google.com/github/minoshim/Num_Analysis/blob/main/FDS/Py/Adv/fds_adv.ipynb
```
と変更すると、Google Colabで開くことができます。<br>
（`https://github.com/`を`https://colab.research.google.com/github/`に変更しています）<br>
また、ノートブックがあるディレクトリのREADMEに設置されている「Open In Colab」ボタンをクリックすることでも開くことができます。

## C++
ローカルのLinux環境（Windows Subsystem for Linuxを含む。多分Macも可）での利用を想定しています。<br>
C++コンパイラとmakeをインストールしてください。<br>
まず、`common`ディレクトリにある`Makefile.inc`で、ご自身の環境に合わせてコンパイラとコンパイラオプションを変更してください。<br>
変更したら、`common`ディレクトリにて
```
cd common
make clean
make
```
とし、`bound.o`が作成されるか確認してください。<br>
ここまでうまくいったら、課題のディレクトリに移動し、
```
cd ../FDS/C++/Adv
make
./a.out
```
とすることで、プログラムが実行されます。上の`FDS/C++/Adv`は例ですので、実行したい課題のディレクトリ（例: `FDS/C++/Diff`, `FVS/C++/Adv`, `Vlasov/Elesta/C++/landau`など）に置き換えてください。<br>
計算結果はPythonで可視化することができます（numpyとmatplotlibが必要）。読み込み用のバッチファイル`batch1d.py`などが準備されていますので、pythonを実行して、
```
exec(open("batch1d.py").read())
```
とします。すると`Input data directory (Ctrl-D to exit):`と聞かれますから、計算を実行したディレクトリであれば`.`と答えて計算結果を読み込みます。2次元Poisson問題では`batch2d.py`、Vlasov問題では`batch_f.py`を使います。<br>
詳しくはバッチファイルの中身をご覧ください。

