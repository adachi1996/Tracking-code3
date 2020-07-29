module const_data
  implicit none
  !=====================================================================================================================共通の定数
  !--------------------------------------------------------------------定数
  double precision,parameter :: pi   =  3.1415926535d0 !円周率
  double precision,parameter :: m0   =  0.5109989461d0 !電子の静止質量エネルギー
  double precision,parameter :: mp0  =  938.2720813d0  !陽子の静止質量エネルギー
  double precision,parameter :: c    =  2.99792458d8 !光速
  double precision,parameter :: L1   =  0.190d0*0.5*0.5!矩形電磁石の幅　(1/2) [m]
  double precision,parameter :: L2   =  0.105d0*0.5       !矩形電磁石の奥行(1/2) [m]
  double precision,parameter :: Ls   =  0.03d0         !FDFの並びのDとFの間の距離 [m]
  !=====================================================================================================================理想磁場軌道計算用
  !------------------------------------------------------------------------------------初期値
  double precision,parameter :: t0   =  0.0d0          !初期時間 [sec]
  double precision,parameter :: z0   =  0.0d0          !初期z [m]
  double precision,parameter :: th0  =  0.0d0          !初期角度 [deg.]
  double precision,parameter :: r0   =  1.0d0          !初期r [m]
  double precision           :: T    =  0.02d0         !初期運動エネルギー [MeV]
  !------------------------------------------------------------------------------------変数θ関係
  double precision,parameter :: dth  =  22.5d0         !保存する角度 [deg.]
  double precision,parameter :: th_h =  1.0d-4         !刻み幅 [deg.]
  !------------------------------------------------------------------------------------ラティス関係
  double precision           :: beF  =  1.74756690474d0*pi/180.0 !2.8125*pi/180.0 !F磁石の見込み角   [deg.]
  double precision           :: beD  =  2.52758977962d0*pi/180.0 !2.8125*pi/180.0 !D磁石の見込み角   [deg.]
  double precision           :: thF  =  22.5d0          *pi/180.0 !F磁石での曲げ角   [deg.]
  double precision           :: M    =  2.0d0              !m値 [1/m]
  double precision,parameter :: N    =  16.0d0             !セル数
  !------------------------------------------------------------------------------------設定値
  double precision,parameter :: eps  =  0.0d0          !漏れ磁場の及ぶ距離 [m]
  !=====================================================================================================================磁場map使用軌道計算用
  !------------------------------------------------------------------------------------初期値
  double precision,parameter :: m_t0   =  0.0d0        !初期時間 [sec]
  double precision           :: m_z0   =  0.0d0        !初期z [m]
  double precision,parameter :: m_th0  =  0.0d0        !初期角度 [deg.]
  double precision,parameter :: m_r0   =  1.0d0!4.60!05843361d0!        !初期r [m]
  double precision           :: m_T    =  0.035d0       !初期運動エネルギー [MeV]
  !------------------------------------------------------------------------------------変数θ関係
  double precision,parameter :: m_dth  =  0.5d0        !保存する角度 [deg.]
  double precision,parameter :: m_th_h =  1.0d-4       !刻み幅 [deg.]
  double precision,parameter :: sym_th =  22.5         !1セルの角度
  !------------------------------------------------------------------------------------閉軌道導出プログラム関係
  double precision,parameter :: cut_z  =  0.305d0        !zの振動の限界値(任意の値)
  !------------------------------------------------------------------------------------読み込む磁場マップ関係
  character(100)             :: file_name1 = 'F[109-095]_171D[085-065]_[FD=273]'
  character(100)             :: file_name2 = '_[11exp(-2z)]_[str+arc]_[type1]'
  character(100)             :: path_name  = 'C:\Users\Kyosuke adachi\Desktop\FORTRAN_SPACE\map_data\'
  character(100)             :: file_name3 = 'Closed_orbit_file_3'
  character(100)             :: path_name2 = 'C:\Users\Kyosuke adachi\Desktop\'
  !=====================================================================================================================共通の設定値
  !------------------------------------------------------------------------------------しらみつぶし
  double precision,parameter :: num  = 10              !回数
  double precision           :: dk   = 0.001d0         !初期刻み幅
  !------------------------------------------------------------------------------------条件
  double precision,parameter :: rev  =  3.0/1.0d0        !周回数
  double precision,parameter :: dr   =  0.0d0          !初期値からのずれ
  double precision,parameter :: dz   =  0.0d0          !初期値からのずれ
  double precision,parameter :: q    = -1.0d0          !電荷(+1 or -1)
  character(100)             :: save_name  = 'test.csv'
  !=====================================================================================================================機能選択
  !トラッキング形式の選択      1:通常トラッキング 2:しらみつぶし       3:アクセプタンス          4:磁場マップ使用    5:m-thF安定領域
  integer,parameter :: select_tracking = 1
  !磁場マップ使用時の機能選択(詳しい説明はmap_trackiを参照)
  !                          1:通常トラッキング 2:複数運動エネルギー  3:磁場マップ生成(座標変換)
  !                          4:垂直FFAG用閉軌道導出計算(位置のみのしらみつぶし)    5:垂直FFAG用閉軌道導出計算(角度含むしらみつぶし)
  !                          6:垂直FFAG用閉軌道導出計算(楕円中心から求めるタイプ)　7:複数エネルギー計算用
  integer,parameter :: select_func     = 1
  !理想的電磁石形状の選択      1:Sector           2:Lectangular        3:Radial sector    4:FDFトリプレット
  integer,parameter :: select_B        = 1
  !FFAGタイプの選択      1:水平FFAG加速器         2:垂直FFAG加速器
  integer,parameter :: select_FFAG     = 2

  !=====================================================================================================================変数
  double precision :: pth                       !角度方向の運動量
  double precision :: r1      ,r2      ,r3      !各軌道位置
  double precision :: beF_rec ,beD_rec ,beS_rec !矩形電磁石の見込み角
  double precision :: roF     ,roD     ,RF      !曲率半径とF磁石での平均半径
  double precision :: BF      ,BD      ,RD      !FとD磁場とD磁石での平均半径
  double precision :: Br      ,Bth     ,Bz      !各磁場成分(r,th,z)
  double precision :: thD     ,thF_rec ,thD_rec !曲げ角
  double precision :: theta   ,r_beta  ,gamma   !1cell中の角度と比速度とγ
  double precision :: r_min   ,th_min  ,z_min   !磁場マップの各データの最小値
  double precision :: drt_r   ,drt_th  ,drt_z   !磁場マップの各データの刻み幅
  double precision :: beta  = (180.0d0/N)*(pi/180.0) !1/2cellの見込み角 [deg.]
  double precision :: temp_th_max !
  double precision :: temp_r_max , temp_z_max , temp_pth_max , temp_pz_max!
  integer          :: check = 0                 !計算終了の判定用
  character        :: filename1*128
  character        :: filename2*128

  !FDFトリプレット用の変数
  double precision :: r0p     ,roDp    ,Ldr     !
  double precision :: Lpp     ,Lppp    !
  double precision :: beF_FDF ,beS_FDF ,beD_FDF !

  !磁場マップ使用時の静止質量エネルギー(水平FFAG加速器は陽子、垂直FFAG加速器は電子)
  double precision :: m_m0

  !計算結果を入れていく配列 [t,th,r,z,pr,pth,pz, Br, Bth, Bz]
  double precision :: particle(10) = (/t0 ,th0 ,r0+dr ,z0+dz ,0.0d0 ,0.0d0 ,0.0d0 ,0.0d0 ,0.0d0 ,0.0d0/)
  double precision :: max_deg = rev*360.0 + th_h !最大角度
  double precision :: h  !刻み幅

  !読み込んだ磁場マップを入れる配列(ポインタ)
  double precision, pointer, dimension(:)       :: r_data
  double precision, pointer, dimension(:)       :: th_data
  double precision, pointer, dimension(:)       :: z_data
  double precision, pointer, dimension(:,:)     :: B_data

end module const_data
