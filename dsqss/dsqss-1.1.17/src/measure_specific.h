//
// measurement info for the SU(2) antiferromagnetic Heisenberg model
// 
// anv  = <Nv>/V
// ene  = ( <Nv> + EBASE ) / V = ( EBASE - T * <Nv> ) / V
// spe  = (<E*E> - <E>^2) / V
// amz  =  <Mz> / V ( should be 0 )
// smz  = <Mz*Mz> / V
// xmz  = <Mz;Mz> / V
// len  = [average worm path]
// xmx  =  <Mx;Mx> / V = [average worm path] * T
// cov  = [coverage of worm paths] = NCYC * [len] / V / BETA
//
// V = Ns = [number of spins]
// EBASE = [sum of base energies of all the interactions]
// Nv = [the number of vertices]
// Mz  = \sum_R (-1)^R Sz(R)
// Mx  = \sum_R (-1)^R Sx(R)
//
// <Q;Q> = \int_0^{\beta} dt <Q(t)Q(0)> / \beta
//       = T * V * \int dt dR <Q(t,R)Q(0,0)>
//

namespace Specific {

  // accumulator specifier
  //    NACC = number of quantities measured at each MC step

  enum                        {  NV1 , EB1 , EB2 , LE1 , MZU1 , MZU2 , MZUA1 , MZUA2 , MZS1 , MZS2 , MZSA1 , MZSA2 , NACC };
  static string ANAME[NACC] = { "nv1","eb1","eb2","le1","mzu1","mzu2","mzua1","mzua2","mzs1","mzs2","mzsa1","mzsa2" };

  // observable specifier
  //    NPHY = number of quantities computed at each set

  enum                        {  ANV , ENE , SPE , LEN , XMX , AMZU , BMZU , SMZU , XMZU , AMZS , BMZS , SMZS , XMZS , NPHY };
  static string PNAME[NPHY] = { "anv","ene","spe","len","xmx","amzu","bmzu","smzu","xmzu","amzs","bmzs","smzs","xmzs"};

}
