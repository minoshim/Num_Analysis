# これは何ですか？
ここは、[計算物理春の学校2026](https://compphysschool.github.io/2026/) 個別講義D-4 「無衝突系の運動論シミュレーション入門」で用いるプログラム置き場です。<br>
- FDS: 有限差分法
- FVS: 有限体積法
- Vlasov: 運動論シミュレーション
- common: 共通C++プログラム

C++とPythonプログラムが準備されていますので、お好きな方をご利用ください。<br>
そのまま使うもよし、機能を追加するもよし、自作の参考にするもよしです。<br>
講義ではPythonを用います。<br>

# どうやって使いますか?
ローカル環境で使う場合は、ご自身のPCにファイルをダウンロードするか、
```
git clone git@github.com:minoshim/Num_Analysis.git
```
してください。<br>
クラウド環境で使う場合は、Google Colabを使いますので、Googleアカウントをご準備ください。
## Python
ディレクトリ`Py`の先にある`main.py`がバッチファイル、`hoge.ipynb`がノートブックです。<br>
バッチファイルをローカル環境で使う場合は、numpy, scipy, matplotlibがインストール済みのpythonを準備します。<br>
pythonを実行して、
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
（`https://github.com/`を`https://colab.research.google.com/github/`に変更しています）

## C++
ローカルのLinux環境（多分Macも可）での利用を想定しています。<br>
c++コンパイラとmakeをインストールしてください。<br>
まず、`common`ディレクトリにある`Makefile.inc`で、ご自身の環境に合わせてコンパイラとコンパイラオプションを変更してください。<br>
変更したら、`common`ディレクトリにて
```
make clean
make
```
とし、`bound.o`が作成されるか確認してください。<br>
ここまでうまくいったら、課題のディレクトリに移動し、
```
make
./a.out
```
とすることで、プログラムが実行されます。<br>
計算結果はPythonで可視化することができます（numpyとmatplotlibが必要）。読み込み用のバッチファイル`batch1d.py`などが準備されていますので、pythonを実行して、
```
exec(open("batch1d.py").read())
```
とします。すると`Input data directory (Ctrl-D to exit):`と聞かれますから、`.`と答えて計算結果を読み込みます。<br>
詳しくはバッチファイルの中身をご覧ください。


