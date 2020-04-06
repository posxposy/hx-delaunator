Haxe port of an incredibly fast [delaunator](https://github.com/mapbox/delaunator) JavaScript library for
[Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) of 2D points.

## Performance Table
Benchmark results of Haxe cross-compilation against original Delaunay JS library (commit `4e6ecd4`)

| &nbsp;           | uniform 100k | gauss 100k | grid 100k | degen 100k | uniform 1&nbsp;million | gauss 1&nbsp;million | grid 1&nbsp;million | degen 1&nbsp;million |
| :--------------- | -----------: | ---------: | --------: | ---------: | ---------------------: | -------------------: | ------------------: | -------------------: |
| **Original Lib** |         79ms |       74ms |      80ms |       33ms |                  1.16s |               1.17ms |              1.02ms |                0.34s |
| Haxe JS          |         69ms |       65ms |      67ms |       37ms |                  1.06s |                1.05s |               0.93s |                0.68s |
| Haxe C++         |         72ms |       73ms |     139ms |       94ms |                  1.10s |                1.06s |               0.40s |                0.62s |
| Haxe C#          |         80ms |       79ms |      71ms |       56ms |                  1.16s |                1.15s |               0.95s |                0.83s |
| Haxe Java        |        118ms |       76ms |      66ms |       42ms |                  1.55s |                1.41s |               1.15s |                0.94s |
| HashLink C       |         94ms |        95s |      86ms |       69ms |                  1.38s |                1.32s |               1.16s |                1.15s |
| HashLink JIT     |        203ms |      197ms |     207ms |      146ms |                  2.63s |                2.74s |               2.52s |                2.63s |

### Performance Comparsion Chart
![100k](https://github.com/dmitryhryppa/hx-delaunator/blob/master/compare_100k.png)
![1m](https://github.com/dmitryhryppa/hx-delaunator/blob/master/compare_1m.png)

### Keep in mind
All these benchmark results depend on the hardware and moon phases.
This comparison was done just for fun and all of these Haxe targets including the original library are incredibly fast enough in the real world.
