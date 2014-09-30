インストール手順
#####################

必要なライブラリ
---------------------

  * MPI
  * BLAS
  * LAPACK


インストール
--------------

.. code-block:: sh

   $ tar xvf dsqss-1.1.16.tar.bz2
   $ cd dsqss-1.1.16
   $ vi runConfigure
        #実行環境にあわせて以下の項目を修正してください。
        #交換法の計算を行うにはモジュールdlaはMPIでのコンパイルが必須です。PARALLEL=YESと指定して下さい
        # CC          :: C compiler command
        # CXX         :: C++ compiler command
        # CXXFLAGS    :: C++ compiler flags
        # CPPFLAGS    :: Preprocessor flags
        # LDFLAGS     :: Linker flags
        # LIBS        :: Libraries(blas,lapack)
        # INSTALL_DIR :: Installation directory
        #* 物性研スパコンシステムBの環境でコンパイルする場合は、
        #* SYSTEM=kashiwaと指定してください。  
   $ ./runConfigure
   $ make
   $ make install

 
$INSTALL_DIR/bin配下に実行モジュール **dla** , **dla_alg** , **lattgene**,  **hamgen_H** が生成されます。

個別モジュールのビルド方法
--------------

.. code-block:: sh

   #lattgene
   $ make lattgene

   #dla_alg
   $ make dla_alg

   #アーカイブの作成　dsqss-x.x.xx.tar.bz2
   $ make dist

   #アーカイブの作成　dsqss-x.x.xx-YYYYmmdd.tar.bz2
   $ make distx

   #アーカイブの作成　dsqss-x.x.xx-YYYYmmddHHMM.tar.bz2
   $ make distx2

   #ユーザマニュアル(html)の生成
   $ make sphinx

   #プログラムドキュメント(html)の生成
   $ make doxygen

アーカイブのバージョン(dsqss-x.x.xx)の変更はconfigure.acのAC_INITを修正してください。


確認テスト
------------

.. code-block:: sh

   $ cd dsqss-1.1.16
   $ ./runConfigure
   $ make
   $ make check

    #テスト後、次のように表示されれば正常にインストールされています
     ++++++++++++++++++++++++++++
     |  dsqss-1.1 Test Passed   |
     ++++++++++++++++++++++++++++

    #インストールが正しくおこなわなければ、次のように表示されます
     ++++++++++++++++++++++++++++
     |  dsqss-1.1 Test Failed   |
     ++++++++++++++++++++++++++++


確認済みの実行環境
------------------


  * PC環境

    * Ubuntu Linux 10.10 64-bit
    * Ubuntu Linux 11.04 64-bit
    * Ubuntu Linux 12.04 64-bit
    * CentOS 5.3 64-bit
     * Gnu C++ Compiler 4.1.x
     * Gnu C++ Compiler 4.4.x
     * Gnu C++ Compiler 4.6.x
     * Intel C++ Compiler 11.1
     * Intel C++ Compiler 12.1
     * OpenMPI 1.4.1
     * OpenMPI 1.4.3
     * Intel MPI 4.0

  * 東京大学物性研究所スパコン

    * `システムB <http://kawashima.issp.u-tokyo.ac.jp/dsqss/sample.html#sample-4-b>`_

  * 京コンピュータ
