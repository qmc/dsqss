/*
 *=================================================================
 *     Header of Class Random for random number generation
 *=================================================================
 *     $Log: Random.h,v $
 *     Revision 1.1  2000/02/10 16:14:16  kenji
 *     Initial revision
 *
 *=================================================================
 *     The copyright holder of the following codes is
 *
 *     Kenji HARADA
 *     Graduate School of Infomatics, Kyoto University,
 *     Kyoto 606-8501, Japan
 *     e-mail: harada@acs.i.kyoto-u.ac.jp
 *     home-page: http://www-fcs.acs.i.kyoto-u.ac.jp/~harada/
 *=================================================================
 */
//! id="$Id: Random.h,v 1.1 2000/02/10 16:14:16 kenji Exp $"
//! author="Kenji Harada"

#ifndef _RANDOM_H_
#include <math.h>
#define IPP 521
#define SIZE_SEED IPP
#define IQQ 32
#define IPQ (IPP-IQQ)
typedef unsigned int Rint;

//: ��������
// ���̃N���X�ł͗������l�n��@��p���Đ������Ă���B
// ���ۂ�t�Ԗڂ̗���X(t)�͎�������v�Z����B<BR>
// X(t) := X(t-32) xor X(t-521)
class Random{
private:
  Rint nrbit,iptr,navr;
  Rint iri[IPP];
  double runit;

private:
  void initialize(Rint irand0,Rint nrbit0);
  //: �����ݒ�
  //!param: irand0 - 521�̏����l�𐶐�����̂ɗp�����
  //!param: nrbit0 - �����̐��x�inrbit0�r�b�g�j
  // �����ݒ������֐�

public:
  Random(Rint *seed,Rint nrbit0);
  //: �����l�Ɛ��x���w�肵�� Class Constructor
  // �����l�z�� seed�Ɋi�[����Ă���521��
  // ���xnrbit0�r�b�g�̏����l�Ƃ��ăZ�b�g����B

  Random(Rint irand0=20000101,Rint nrbit0=32);
  //: ��A���x�w�肵�� Class Constructor
  // �� irand0���g����nrbit0�r�b�g�̐��x���������̂��߂̏����ݒ������B

  void setSeed(Rint irand0,Rint nrbit0=32);
  //: ��A���x���w�肵���ď����ݒ�
  //!param: irand0 - 521�̏����l�𐶐�����̂ɗp�����
  //!param: nrbit0 - �����̐��x�inrbit0�r�b�g�j
  // ���irand0�ɐ��xnrbit0�r�b�g�ōď����ݒ������B

  void setSeed(Rint *seed,Rint nrbit0);
  //: �����l�Ɛ��x���Đݒ�
  //!param: seed - 521�̏����l���i�[����Ă���L���̈�ւ̃|�C���^
  //!param: nrbit0 - �����̐��x�inrbit0�r�b�g�j
  // �����l�z�� init�Ɋi�[����Ă���521��
  // ���xnrbit0�r�b�g�̏����l�Ƃ��ăZ�b�g����B

  Rint getSeed (Rint *seed);
  //: ���݂̗����̏����l�Ɛ��x��Ԃ�
  //!param: seed - ���݂̗����̏����l(521��)���i�[����L���̈�ւ̃|�C���^
  // ���݂̗����̏����l�i521�j���|�C���^seed�̎w���L���̈�Ɋi�[����B
  // ����ɐ��x�̃r�b�g����߂�l�Ƃ��ĕԂ��B

  double Uniform(void);

  Rint Int(Rint);
  //: ���U��l���z (0, 1, .. ,ilimit-1)
  //!param: ilimit - �͈͂̏��
  // �߂�l�� 0�ȏ�Ailimit�����̗���(����)��Ԃ��B

  Rint Int(void);

  void Uniform(Rint nr,Rint *ir);
  //: ���U��l���z����(�����l�z��)
  //!param: nr - �������闐���̌�
  //!param: ir - �������������l�����i0�ȏ�A2**nrbit0-1 �����j���i�[����Ԓn�ւ̃|�C���^
  // �����l�����i0�ȏ�A2**nrbit0-1�����j��nr�������A�����|�C���^ ir�̎w���Ԓn���珇�ԂɊi�[����B

  void Uniform(Rint nr,double *rx);
  //: ��l���z����(�����l�z��)
  //!param: nr - �������闐���̌�
  //!param: rx - �������������l�����i0�ȏ�A1�����j���i�[����Ԓn�ւ̃|�C���^
  // �����l�����i0�ȏ�A1�����j��nr�������A�����|�C���^ rx�̎w���Ԓn���珇�ԂɊi�[����B

  void Int(Rint nr,Rint *ir,Rint ilimit);
  //: ����t�����U��l���z(�����l�z��)
  //!param: nr - �������闐���̌�
  //!param: ir - �������������l�����i0�ȏ�Ailimit �����j���i�[����Ԓn�ւ̃|�C���^
  //!param: ilimit - �������鐮���l�����̏��
  // �����l�����i0�ȏ�Ailimit�����j��nr�������A�����|�C���^ ir�̎w���Ԓn���珇�ԂɊi�[����B

  double Gauss(){
    double theta;
    theta = 6.283185307179586477*Uniform();
    return sqrt(-2e0*log(1e0-Uniform()))*sin(theta);
  }
  //: ���K���z

  double Exp(){return -log(1e0-Uniform());}
  //: ���K�w�����z

  int Binary(double P){return ((int) (Uniform()/P));}
  //: �Q�����z

  void Perm(Rint, int*);

  void Scramble(Rint, int*);

};

#define _RANDOM_H_
#endif
