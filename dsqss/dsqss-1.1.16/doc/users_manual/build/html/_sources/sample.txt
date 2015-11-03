計算例
#########

sample-1 : S=1/2 ハイゼンベルグモデル
--------------------------------------------
もっとも簡単な例として、S=1/2 のスピン２個だけからなる系を考えます。 以下は、これらが反強磁性的に結合していて、磁場がかかっていない場合の計算例です。 アルゴリズムファイル（ワームの散乱確率や、バーテックス密度を含む）や 格子ファイル（単位胞の情報を含む）を直接編集して、dla 単体で 実行することも可能ですが、ここでは、hamgen_H, dla_alg, latgene を利用して、 これらのファイルを自動生成して計算する例を示します。 dla.exe, hamgen_H.exe, dla_alg.exe, latgene.exe などの 実行ファイルのあるディレクトリにはパスが通っているとします。 （通っていなければ、パスも併せて指定してください。） 各実行ファイルの関係や入出力ファイルの詳細については、 実行方法ページ を参照。

入力ファイル（sample.inp）例:
 .. code-block:: sh

    runtype = 0                # シリアル実行

    nmcse = 1000               # 平衡化のために捨てられる初期モンテカルロステップ数
    nmcsd = 1000               # 次setとのあいだの相関を減らすために捨てられる set 間インターバルのステップ数
    nmcs  = 10000              # 各setごとの物理量測定ステップ数
    nset  = 10                 # set数
    seed  = 31415              # 乱数シード
    nvermax = 10000            # バーテックス数上限（実行中に上限に達すると自動終了）
    nsegmax = 10000            # セグメント数上限（実行中に上限に達すると自動終了）

    algfile = algorithm.xml    # アルゴリズムファイル名
    latfile = lattice.xml      # 格子ファイル名
    outfile = sample.log       # 出力ファイル名

計算実行例 :
 .. code-block:: sh

    $ $INSTALL_DIR/hamgen_H 1 -0.5 0.0          # ハミルトニアン生成．パラメータはスピン長（2*S），交換相互作用（J），外部磁場（F）の順
                                                # 長さ２の周期境界１次元系として実行するので，実際の交換相互作用は２倍の -1.0 である．
                                                # また，F は全ハミルトニアンをペアハミルトニアンに分割したときの個々のペアハミルトニアン
                                                # の中での磁場なので，d 次元立方格子の場合は，「通常の」１スピンあたりの磁場は 2dF になる．
    $ $INSTALL_DIR/dla_alg                      # ハミルトニアンから，アルゴリズム生成．
    $ $INSTALL_DIR/lattgene 1 2 1.0             # 格子生成．パラメータは次元（D），サイズ（L），逆温度（\beta）の順．
    $ mpirun -np 1 $INSTALL_DIR/dla sample.inp  # 並列度１で実行．

　サンプルスクリプトの使用例 :
 .. code-block:: sh

  $ $INSTALL_DIR/samp/run_test.sh
    # mpi,実行モジュールのパス修正が必要です。
    定義生成ファイルの作成から本計算まで対話形式で実行できます。

出力ファイル（sample.log.000）例:
 .. code-block:: sh

    C This is DSQSS ver.1.1

    P D       =            1                         # 次元
    P L       =            2                         # １辺の長さ
    P BETA    =       1.0000000000000000             # 逆温度
    P NSET    =           10                         # セット数
    P NMCSE   =         1000                         # 平衡化ステップ数
    P NMCSD   =         1000                         # セット間インターバル
    P NMCS    =        10000                         # 測定ステップ数
    P SEED    =        31415                         # 乱数シード
    P NSEGMAX =        10000                         # 最大セグメント数
    P NVERMAX =        10000                         # 最大バーテックス数
    P NCYC    =            2                         # １ステップあたりのワームサイクル数（プログラムが決定）
    P ALGFILE = algorithm.xml                        # アルゴリズムファイル名
    P LATFILE = lattice.xml                          # 格子ファイル名
    P OUTFILE = sample.log.000                       # 出力ファイル名
    R anv    = 7.5180000000e-02 7.2092533132e-04     # 物理量（R "ラベル" = "平均値" "統計誤差"） 
    R ene    = -1.1287250000e-01 1.1319872521e-03
    R spe    = 1.2403653237e-01 1.1809645533e-03
    R len    = 1.2023753000e+00 6.6838265039e-04
    R xmx    = 3.0059382501e-01 1.6709566260e-04
    R amzu   = 4.9500000000e-04 1.0936800771e-03
    R bmzu   = 4.9500000000e-04 1.0936800771e-03
    R smzu   = 1.7461500000e-01 1.0417413307e-03
    R xmzu   = 1.7461500000e-01 1.0417413307e-03
    R amzs   = 1.6500000000e-04 1.0586167599e-03
    R bmzs   = -1.9144707385e-04 9.9337993655e-04
    R smzs   = 3.2538500000e-01 1.0417413307e-03
    R xmzs   = 3.0049234531e-01 9.0849325770e-04
    I [the maximum number of segments]          = 17    # 実際の最大セグメント数 
    I [the maximum number of vertices]          = 10    # 実際の最大バーテックス数
    I [the maximum number of reg. vertex info.] = 3  

sample-2 : レプリカ交換法を用いた拡張アンサンブル計算(磁場)
------------------------------------------------------------

DSQSSでは、磁場と逆温度のパラメータを対象としレプリカ交換法を用いた拡張アンサンブル計算が可能です。
sample-2では磁場をパラメータとしたレプリカ交換法を用いた拡張アンサンブル計算の設定、計算方法を示します。レプリカ数8、磁場は0.4から0.02間隔、最大交換数100の系を例にします。
各レプリカ[0-7]において、[0.4,0.42,0.44,0.44・・・0.54]が磁場の初期値としてセットされます。

計算に先立って、アルゴリズム定義ファイル、格子定義ファイルを生成します。各レプリカによって磁場の初期値が異なるため、磁場の値に対応するアルゴリズム定義ファイルを生成する必要があります。各磁場値にて生成されたアルゴリズム定義ファイルalgorithm.xmlをレプリカ番号.xmlの名前で保存してください。また、レプリカ数分のアルゴリズム定義ファイルの他にダミーのアルゴリズムファイル定義が必要になります。入力ファイルの変数algfile(デフォルトはalgorthm.xml)で指定されたファイルを用意してください。

　定義ファイル生成例 :
 .. code-block:: sh

  $ $INSTALL_DIR/bin/lattgene 1 4 0.1
    # === レプリカ０ ===
  $ $INSTALL_DIR/bin/hamgen_H 1 1.0 0.4
  $ $INSTALL_DIR/bin/dla_alg 
  $ mv algorithm.xml 0.xml

    # === レプリカ１ ===
  $ $INSTALL_DIR/bin/hamgen_H 1 1.0 0.42
  $ $INSTALL_DIR/bin/dla_alg 
  $ mv algorithm.xml 1.xml
       .
       .
    # === レプリカ７ ===
  $ $INSTALL_DIR/bin/exact_H 1 1.0 0.54
  $ $INSTALL_DIR/bin/dla_alg
  $ mv algorithm.xml 7.xml

    # === ダミーファイルの作成 ===
  $ cp 0.xml algorithm.xml(0.xmlは任意)

次に、計算制御パラメータの ``runtype`` , ``nset`` , ``nrep`` , ``vf`` , ``df`` を修正します。

　入力ファイル例 :
 .. code-block:: bash

      # ==== 入力ファイル例  qmc.inp  ====
      # == RUNTYPE ==
      runtype = 1 #磁場によるレプリカ交換計算

      # == PARAMETER ==
      nmcse   = 10
      nmcsd   = 100
      nmcs    = 500
      nset    = 100    #最大交換数(交換判定回数)
      seed    = 71314416
      nrep    =   8    #レプリカ数
      vf      = 0.4    # 磁場の最小値
      df      = 0.02   # 各レプリカの磁場の間隔
      nvermax = 10000
      nsegmax = 10000

      # == OUTPUT_FILE ==
      outfile = qmc.log

レプリカ数＝MPI並列数として、計算を実行してください。

　計算実行例 :
 .. code-block:: sh

  $ mpirun -np 8 $INSTALL_DIR/dla qmc.inp
     # qmc.log.000-007が生成されます。

　サンプルスクリプトの使用例 :
 .. code-block:: sh

  $ $INSTALL_DIR/samp/run_test.sh
    # mpi,実行モジュールのパス修正が必要です。
    定義生成ファイルの作成から本計算まで対話形式で実行できます。


　計算結果例 : qmc.log.002
 .. code-block:: sh

    C This is DSQSS ver.1.1

    P D       =            1
    P L       =            4
    P BETA    =       0.1000000000000000
    P NSET    =          100
    P NMCSE   =           10
    P NMCSD   =          100
    P NMCS    =          500
    P SEED    =     71314416
    P NSEGMAX =        10000
    P NVERMAX =        10000
    P NCYC    =            4
    P ALGFILE = algorithm.xml
    P LATFILE = lattice.xml
    P OUTFILE = qmc.log.002
    R anv    = 1.1500000000e-03 1.1135075168e-04
    R ene    = -1.2234756810e-01 2.7178997390e-03
    R spe    = 1.3032337345e-02 1.2234719739e-04
    R len    = 1.0466299052e-01 7.1219677615e-05
    R xmx    = 2.6165747630e-01 1.7804919404e-04
    R amzu   = 5.1245000000e-02 1.2152959277e-03
    R bmzu   = 5.1245000000e-02 1.2152959277e-03
    R smzu   = 2.6935500000e-01 1.4698707254e-03
    R xmzu   = 2.6935500000e-01 1.4698707254e-03
    R amzs   = -2.5350000000e-03 1.1876097504e-03
    R bmzs   = -2.6581629401e-03 1.1902119933e-03
    R smzs   = 2.3599500000e-01 1.3198847325e-03
    R xmzs   = 2.3551502705e-01 1.3173450648e-03
    I [the maximum number of segments]          = 13
    I [the maximum number of vertices]          = 9
    I [the maximum number of reg. vertex info.] = 3
    I [the index of replica exchange ]         = 5

計算結果ファイルqmc.log.002は、プロセサー2が最終的に担当するレプリカの計算結果を出力します。
ファイル最後のI [the index of replica exchange ] の行にある数値（ここでは5）が、最初に設定したレプリカの認識番号です。すなわち、この出力ファイルはレプリカ番号5（磁場0.5)の計算結果になります。概念図を次にしめします。


.. image :: replica_fig2.png



sample-3 :  レプリカ交換法を用いた拡張アンサンブル計算(逆温度)
--------------------------------------------------------------------

sample-3では逆温度をパラメータとしたレプリカ交換法を用いた拡張アンサンブル計算の設定、計算方法を示します。レプリカ数8、逆温度は0.1から0.12間隔、最大交換数100の系を例にします。
各レプリカ[0-7]において、[0.1,0.22,0.34,0.46・・・0.94]が逆温度の初期値としてセットされます。

計算に先立って、アルゴリズム定義ファイル、格子定義ファイルを生成します。
詳しくは、「実行手順」の「実行方法」を参照してください。



　定義ファイル生成例 :
 .. code-block:: sh

  $ $INSTALL_DIR/hamgen_H 1 1.0 0.4
  $ $INSTALL_DIR/dla_alg 
       # algorithm.xmlが生成されます。
  $ $INSTALL_DIR/lattgene 1 4 0.1(注)
       # (注)Tの値は計算制御ファイルで与えるので、ここでは任意の値で構いません。
       # lattice.xmlが生成されます。

各レプリカに対応するアルゴリズム定義ファイルを生成する必要があります。各レプリカでの磁場値は等しくアルゴリズム定義ファイルは同じですが、アルゴリズム定義ファイルalgorithm.xmlをレプリカ番号.xmlの名前で保存してください。また、レプリカ数分のアルゴリズム定義ファイルの他にダミーのアルゴリズムファイル定義が必要になります。入力ファイルの変数algfile(デフォルトはalgorthm.xml)で指定されたファイルを用意してください。

次に、計算制御パラメータの ``runtype`` , ``nset`` , ``nrep`` , ``vb`` , ``db`` を修正します。

　入力ファイル例 :
 .. code-block:: bash


      # ==== 入力ファイル例  qmc.inp  ====
      # == RUNTYPE ==
      runtype = 2 #逆温度によるレプリカ交換計算

      # == PARAMETER ==
      nmcse   = 10
      nmcsd   = 100
      nmcs    = 500
      nset    = 100    #最大交換数(交換判定回数)
      seed    = 71314416
      nrep    =   8    #レプリカ数
      vb      = 0.1    # 逆温度の最小値
      db      = 0.1    # 各レプリカの逆温度の間隔
      nvermax = 10000
      nsegmax = 10000

      # == OUTPUT_FILE ==
      outfile = qmc.log

レプリカ数＝MPI並列数として、計算を実行してください。

　計算実行例 :
 .. code-block:: sh

  $ mpirun -np 8 $INSTALL_DIR/dla qmc.inp
     # qmc.log.000-007が生成されます。


　サンプルスクリプトの使用例 :
 .. code-block:: sh

  $ $INSTALL_DIR/samp/run_test.sh 
    # mpi,実行モジュールのパス修正が必要です。
    定義生成ファイルの作成から本計算まで対話形式で実行できます。

　計算結果例 : qmc.log.006
 .. code-block:: sh

    C This is DSQSS ver.1.1
    
    P D       =            1
    P L       =            4
    P BETA    =       0.3400000000000000
    P NSET    =          100
    P NMCSE   =           10
    P NMCSD   =          100
    P NMCS    =          500
    P SEED    =     71314416
    P NSEGMAX =        10000
    P NVERMAX =        10000
    P NCYC    =            4
    P ALGFILE = algorithm.xml
    P LATFILE = lattice.xml
    P OUTFILE = qmc.log.006
    R anv    = 1.2690000000e-02 3.5009233991e-04
    R ene    = -3.2201857198e-01 2.5566094999e-03
    R spe    = 1.1614163054e-01 4.9001930595e-04
    R len    = 3.8298037722e-01 3.8403676477e-04
    R xmx    = 2.8160321854e-01 2.8237997409e-04
    R amzu   = 1.5311500000e-01 1.2735033799e-03
    R bmzu   = 1.5311500000e-01 1.2735033799e-03
    R smzu   = 3.5586500000e-01 1.7006364673e-03
    R xmzu   = 3.5586500000e-01 1.7006364673e-03
    R amzs   = -1.0150000000e-03 1.0338330903e-03
    R bmzs   = -9.9371314432e-04 1.0226770536e-03
    R smzs   = 1.9676500000e-01 1.1740342044e-03
    R xmzs   = 1.9279519342e-01 1.1575030817e-03
    I [the maximum number of segments]          = 21
    I [the maximum number of vertices]          = 13
    I [the maximum number of reg. vertex info.] = 3
    I [the index of replica exchange ]         = 2

計算結果ファイルqmc.log.006は、プロセサー6が最終的に担当するレプリカの計算結果を出力します。
ファイル最後のI [the index of replica exchange ] の行にある数値（ここでは2）が、最初に設定したレプリカの認識番号です。すなわち、この出力ファイルはレプリカ番号2（逆温度0.34)の計算結果になります。概念図を次にしめします。



.. image :: replica_fig1.png

sample-4 : 物性研システムBでの利用方法（インタラクティブノード）
----------------------------------------------------------------
DSQSSは物性研スパコンにあらかじめインストールされています。
下記のようにパスを通して使用してください。

 .. code-block:: sh

  $ export PATH=/opt/nano/DSQSS/dsqss-1.1.15/bin

アルゴリズムや格子の定義方法は通常の方法と同じです。前述の実行サンプルを参照してください。

　モンテカルロ計算dlaの実行例 :
 .. code-block:: sh

  $ bsub -q i32 -n 1 -W 1 mpijob $DSQSS_PATH/dla qmc.inp>qmc.log
  # 並列数4で実行。詳細情報はシステムBの利用マニュアルを参照してください。

sample-5 : 物性研システムBでの利用方法（バッチジョブとして実行する場合）
------------------------------------------------------------------------
システムBでキューを用いてバッチ処理をおこなう場合は、次のようなbsubスクリプトを作成し、bsubに投入してください。

　run.sh :B1キューを用いて4並列計算をおこなう場合の例
 .. code-block:: sh

  #!/bin/bash

  #bsub -q B1
  #bsub -n 1
  #bsub -w 1
 
  #bsub -J test
  #bsub -o test.log
  #bsub -N
  #bsub -B

  mpijob /opt/nano/DSQSS/dsqss-1.1.15/bin/dla qmc.inp


 .. code-block:: sh

  $ bsub < run.sh
