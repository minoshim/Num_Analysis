# これは何ですか？
ここは、[計算物理春の学校2026](https://compphysschool.github.io/2026/) 個別講義d-4 「無衝突系の運動論シミュレーション入門」で用いるプログラム置き場です。<br>
- FDS: 有限差分法
- FVS: 有限体積法
- Vlasov: 運動論シミュレーション
- common: 共通C++プログラム

C++とPythonプログラムが配布されていますので、お好きな方をお使いください。<br>
そのまま使うもよし、自作の参考にするもよしです。<br>
講義ではPythonを用います。<br>

# どうやって使いますか?
ローカル環境で使う場合は、ファイルをダウンロードするか、`> git clone git@github.com:minoshim/Num_Analysis.git`してください。<br>
クラウド環境で使う場合は、Google Colabで開きますので、Googleアカウントをご準備ください。
## Python
フォルダ`Py`の先にある`main.py`がバッチファイル、`hoge.ipynb`がノートブックです。<br>
バッチファイルを使う場合は、numpy, scipy, matplitlibがインストール済みのpython環境を使います。
```
> python
>>>  exec(open("main.py").read())
```
で実行されます。<br>
詳しくはバッチファイルの中身をご覧ください。

ノートブックを使う場合は、例えば
```
https://github.com/minoshim/Num_Analysis/blob/main/FDS/Py/Adv/fds_adv.ipynb
```
のURLを
```
https://colab.research.google.com/github/minoshim/Num_Analysis/blob/main/FDS/Py/Adv/fds_adv.ipynb
```
と変更すると、Google Colabで開くことができます。<br>
（`https://github.com/`を`https://colab.research.google.com/`に変更しています）

## C++
Linux環境（多分Macも可）での利用を想定しています。<br>
c++コンパイラとmakeをインストールしてください。<br>
まず、`common`ディレクトリにある`Makefile.inc`で、ご自身の環境に合わせてコンパイラとコンパイラオプションを変更してください。<br>
変更したら、`common`ディレクトリにて
```
> make clean
> make
```
とし、`bound.o`が作成されるか確認してください。<br>
ここまでうまくいったら、課題のディレクトリまで移動し、
```
> make
> ./a.out
```
とすることで、プログラムが実行されます。<br>
実行結果はPythonで可視化することができます。読み込み用のバッチファイル`batch1d.py`などが準備されていますので、
```
> python
>>> exec(open("batch1d.py").read())
```
とします。すると`Input data directory (Ctrl-D to exit):`と聞かれますから、`.`と答えることで計算結果を読み込みます。<br>
詳しくはバッチファイルの中身をご覧ください。


