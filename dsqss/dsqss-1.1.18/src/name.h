
#ifndef NAME_H
#define NAME_H

namespace STYPE { int UNDEF = -1; }
namespace HTYPE { int UNDEF = -1; int ABSENT = -2; }
namespace TTYPE { int UNDEF = -1; }
namespace ITYPE { int UNDEF = -1; }
namespace VTYPE { int UNDEF = -1; }
namespace DIR   { int UNDEF = -1; }
namespace STATE { int UNDEF = -1; }
namespace ICONF { int UNDEF = -1; }
namespace ICC   { int UNDEF = -1; }
namespace UORD  { int UNDEF = -1; int UP = 0; int DOWN = 1; }

namespace VCAT  { 
  int UNDEF = -1; 
  int TERM = 0; 
  int WORM = 1; 
  int INT  = 2; 
}

#define INF 1.0e+14
#define EPS 1.0e-14

#endif
