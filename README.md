# benchmark_calc_check_degit

C++でマイナンバーのチェックデジットを計算する

- [C++でマイナンバーのチェックデジットを計算する - Qiita](http://qiita.com/yumetodo/items/600ca0df422010cbc4c1)
- [SIMD intrinsicでチェックディジットを計算してみる - Qiita](http://qiita.com/YSRKEN/items/4ca7229c98640a71bdad)
- [Calculating Japan's My Number check digit using AVX intrinsics](https://gist.github.com/MaverickTse/b78eff8fcc70962e0ee7a21b985bbaa9)

を検証・ベンチマークするコード群です。

## Dependency

- [Sprout](https://github.com/bolero-MURAKAMI/Sprout)

```
git clone --recursive https://github.com/yumetodo/benchmark_calc_check_degit.git
```

or

```
git clone https://github.com/yumetodo/benchmark_calc_check_degit.git
git submodule update --init
```
## Supported Compiler

- Visual Studio 2015 Update3
- Visual Studio 2017
