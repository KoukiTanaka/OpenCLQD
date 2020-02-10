# OpenCLQD
OpenCL上でQDライブラリを利用するためのプログラムです。

## Description
このプログラムは [David H. Bailey](https://www.davidhbailey.com/dhbsoftware/)らが実装した4倍精度浮動小数点数(dd_real)と8倍精度浮動小数点数(qd_real)を扱うネイティブC++用のライブラリであるQDライブラリをOpenCLのカーネルコード内で利用することができるように、OpenCL C言語に移植したプログラムです。

## Usage
カーネルコード内でincludeをするだけで利用することができます。

また、ホスト側とデバイス側での型の関係はdd_real => cl_dd , qd_real => cl_qd　となっています。また、特別な操作を行うことなくOpenCL APIの関数でデータのやり取りを行うことができます。
