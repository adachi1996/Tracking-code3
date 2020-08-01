module const_data
  implicit none
  !=====================================================================================================================���ʂ̒萔
  !--------------------------------------------------------------------�萔
  double precision,parameter :: pi   =  3.1415926535d0 !�~����
  double precision,parameter :: m0   =  0.5109989461d0 !�d�q�̐Î~���ʃG�l���M�[
  double precision,parameter :: mp0  =  938.2720813d0  !�z�q�̐Î~���ʃG�l���M�[
  double precision,parameter :: c    =  2.99792458d8 !����
  double precision,parameter :: L1   =  0.190d0*0.5*0.5!��`�d���΂̕��@(1/2) [m]
  double precision,parameter :: L2   =  0.105d0*0.5       !��`�d���΂̉��s(1/2) [m]
  double precision,parameter :: Ls   =  0.03d0         !FDF�̕��т�D��F�̊Ԃ̋��� [m]
  !=====================================================================================================================���z����O���v�Z�p
  !------------------------------------------------------------------------------------�����l
  double precision,parameter :: t0   =  0.0d0          !�������� [sec]
  double precision,parameter :: z0   =  0.0d0          !����z [m]
  double precision,parameter :: th0  =  0.0d0          !�����p�x [deg.]
  double precision,parameter :: r0   =  1.0d0          !����r [m]
  double precision           :: T    =  0.02d0         !�����^���G�l���M�[ [MeV]
  !------------------------------------------------------------------------------------�ϐ��Ɗ֌W
  double precision,parameter :: dth  =  22.5d0         !�ۑ�����p�x [deg.]
  double precision,parameter :: th_h =  1.0d-4         !���ݕ� [deg.]
  !------------------------------------------------------------------------------------���e�B�X�֌W
  double precision           :: beF  =  1.74756690474d0*pi/180.0 !2.8125*pi/180.0 !F���΂̌����݊p   [deg.]
  double precision           :: beD  =  2.52758977962d0*pi/180.0 !2.8125*pi/180.0 !D���΂̌����݊p   [deg.]
  double precision           :: thF  =  22.5d0          *pi/180.0 !F���΂ł̋Ȃ��p
  double precision           :: M    =  2.0d0              !m�l [1/m]
  double precision,parameter :: N    =  16.0d0             !�Z����
  !------------------------------------------------------------------------------------�ݒ�l
  double precision,parameter :: eps  =  0.0d0          !�R�ꎥ��̋y�ԋ��� [m]
  !=====================================================================================================================����map�g�p�O���v�Z�p
  !------------------------------------------------------------------------------------�����l
  double precision,parameter :: m_t0   =  0.0d0        !�������� [sec]
  double precision,parameter :: m_th0  =  0.0d0        !�����p�x [deg.]
  double precision,parameter :: m_r0   =  0.99756891d0!4.60!05843361d0!        !����r [m]
  double precision           :: m_z0   =  0.14116665d0        !����z [m]
  double precision           :: m_T    =  0.026d0       !�����^���G�l���M�[ [MeV]
  integer         ,parameter :: select_initial = 1     !���ˊp�x���w�肷��ꍇ��1�ɂ��Ĉȉ�(Pr, Pth, Py)����͂���B��������Ȃ��Ȃ�0
                                                       !�����Fm=2, alpha=1.71, T=26 [keV]�̓��ˏ���
  double precision,parameter :: m_pr0  =  0.44934438d-4        !����Pr  [MeV/c]
  double precision,parameter :: m_pth0 =  0.16506950d0        !����Pth [MeV/c]
  double precision,parameter :: m_py0  =  0.55750822d-4        !����Py  [MeV/c]
  !------------------------------------------------------------------------------------�ϐ��Ɗ֌W
  double precision,parameter :: m_dth  =  22.5d0        !�ۑ�����p�x [deg.]
  double precision,parameter :: m_th_h =  1.0d-4       !���ݕ� [deg.]
  double precision,parameter :: sym_th =  22.5         !1�Z���̊p�x
  !------------------------------------------------------------------------------------�����֌W
  double precision,parameter :: th_RF  =  180.0d0+5.625d0         !�����󓴂̐ݒu�p�x(��)
  double precision,parameter :: th_s   =  16.875d0*pi/180.0         !�i�s�����֍��W�ϊ����邽�߂̊p�x(��s)
  double precision,parameter :: RF_kV  =  0.0d0         !�����d�� [kV]
  !------------------------------------------------------------------------------------�O�����o�v���O�����֌W
  double precision,parameter :: cut_z  =  0.24d0        !z�̐U���̌��E�l(�C�ӂ̒l)
  !------------------------------------------------------------------------------------�ǂݍ��ގ���}�b�v�֌W
  character(100)             :: file_name1 = 'F[109-095]_175D[085-065]_[FD=273]'
  character(100)             :: file_name2 = '_[11exp(-2z)]_[str+arc]_[type1]'
  character(100)             :: path_name  = 'C:\Users\Kyosuke adachi\Desktop\FORTRAN_SPACE\map_data\'
  character(100)             :: file_name3 = 'Closed_orbit_file_3'
  character(100)             :: path_name2 = 'C:\Users\Kyosuke adachi\Desktop\'
  !=====================================================================================================================���ʂ̐ݒ�l
  !------------------------------------------------------------------------------------�ȉ~���S���玟�̓��ˏ��������߂��
  integer,parameter          :: set_num  = 10            !��
  double precision,parameter :: set_T0   = 20.0d0
  double precision,parameter :: set_T1   = 70.0d0
  double precision,parameter :: set_dT   = 2.0d0
  double precision,parameter :: set_rev  = 50.0d0
  character(100)             :: set_name = '20200801_m=2'
  !------------------------------------------------------------------------------------����݂Ԃ�(�g���Ă��Ȃ�)
  double precision           :: dk   = 0.001d0         !�������ݕ�
  !------------------------------------------------------------------------------------����
  double precision,parameter :: rev  =  2.0/1.0d0      !����
  double precision,parameter :: dr   =  0.0d0          !�����l����̂���
  double precision,parameter :: dz   =  0.0d0          !�����l����̂���
  double precision,parameter :: q    = -1.0d0          !�d��(+1 or -1)
  character(100)             :: save_name  = 'test.csv'
  !=====================================================================================================================�@�\�I��
  !�g���b�L���O�`���̑I��      1:�ʏ�g���b�L���O 2:����݂Ԃ�       3:�A�N�Z�v�^���X          4:����}�b�v�g�p    5:m-thF����̈�
  integer,parameter :: select_tracking = 4
  !����}�b�v�g�p���̋@�\�I��(�ڂ���������map_tracki���Q��)
  !                          1:�ʏ�g���b�L���O 2:�����^���G�l���M�[  3:����}�b�v����(���W�ϊ�)
  !                          4:����FFAG�p�O�����o�v�Z(�ʒu�݂̂̂���݂Ԃ�)    5:����FFAG�p�O�����o�v�Z(�p�x�܂ނ���݂Ԃ�)
  !                          6:����FFAG�p�O�����o�v�Z(�ȉ~���S���狁�߂�^�C�v)�@7:�����G�l���M�[�v�Z�p
  integer,parameter :: select_func     = 6
  !���z�I�d���Ό`��̑I��      1:Sector           2:Lectangular        3:Radial sector    4:FDF�g���v���b�g
  integer,parameter :: select_B        = 1
  !FFAG�^�C�v�̑I��      1:����FFAG������         2:����FFAG������
  integer,parameter :: select_FFAG     = 2

  !=====================================================================================================================�ϐ�
  double precision :: pth                       !�p�x�����̉^����
  double precision :: r1      ,r2      ,r3      !�e�O���ʒu
  double precision :: beF_rec ,beD_rec ,beS_rec !��`�d���΂̌����݊p
  double precision :: roF     ,roD     ,RF      !�ȗ����a��F���΂ł̕��ϔ��a
  double precision :: BF      ,BD      ,RD      !F��D�����D���΂ł̕��ϔ��a
  double precision :: Br      ,Bth     ,Bz      !�e���ꐬ��(r,th,z)
  double precision :: thD     ,thF_rec ,thD_rec !�Ȃ��p
  double precision :: theta   ,r_beta  ,gamma   !1cell���̊p�x�Ɣ䑬�x�ƃ�
  double precision :: r_min   ,th_min  ,z_min   !����}�b�v�̊e�f�[�^�̍ŏ��l
  double precision :: drt_r   ,drt_th  ,drt_z   !����}�b�v�̊e�f�[�^�̍��ݕ�
  double precision :: beta  = (180.0d0/N)*(pi/180.0) !1/2cell�̌����݊p [deg.]
  double precision :: temp_th_max !
  double precision :: temp_r_max , temp_z_max , temp_pth_max , temp_pz_max!
  integer          :: check = 0                 !�v�Z�I���̔���p
  character        :: filename1*128
  character        :: filename2*128

  !FDF�g���v���b�g�p�̕ϐ�
  double precision :: r0p     ,roDp    ,Ldr     !
  double precision :: Lpp     ,Lppp    !
  double precision :: beF_FDF ,beS_FDF ,beD_FDF !

  !����}�b�v�g�p���̐Î~���ʃG�l���M�[(����FFAG������͗z�q�A����FFAG������͓d�q)
  double precision :: m_m0

  !�v�Z���ʂ����Ă����z�� [t,th,r,z,pr,pth,pz, Br, Bth, Bz]
  double precision :: particle(10) = (/t0 ,th0 ,r0+dr ,z0+dz ,0.0d0 ,0.0d0 ,0.0d0 ,0.0d0 ,0.0d0 ,0.0d0/)
  double precision :: max_deg = rev*360.0d0 !�ő�p�x
  double precision :: h  !���ݕ�

  !�ǂݍ��񂾎���}�b�v������z��(�|�C���^)
  double precision, pointer, dimension(:)       :: r_data
  double precision, pointer, dimension(:)       :: th_data
  double precision, pointer, dimension(:)       :: z_data
  double precision, pointer, dimension(:,:)     :: B_data

end module const_data
