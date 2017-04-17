
// Known bugs:
// (1) Tags must be separated from other words by one or more spaces or 
//     line breaks. Otherwise they are not recognized.
// (2) The comment identifiers <!-- and --> must be separated by spaces or
//     line breaks, too. 

#ifndef XML_H
#define XML_H

#include "io.h"
#include "array.h"

namespace XML {

//#######################################################################

class Block;

//======================================================================

inline bool isCommentOpening(const string& key) {
  string sopen = "<!--";
  int n = key.length();
  if ( n < 4 ) return false;
  if ( key.substr(0,4) == sopen ) return true;
  return false;
}

//======================================================================

inline bool isOpening(const string& key) {
  int n = key.length();
  if ( key[0] != '<' ) return false;
  if ( key[1] == '/' ) return false;
  if ( key[n-1] != '>' ) return false;
  return true;
}

//======================================================================

inline bool isClosing(const string& key) {
  int n = key.length();
  if ( key[0] != '<' ) return false;
  if ( key[1] != '/' ) return false;
  if ( key[n-1] != '>' ) return false;
  return true;
}

//======================================================================

inline const string getOpeningName(const string& key) {
  if ( !isOpening(key) ) exit(0);
  int n = key.length();
  string ans = key.substr(1,n-2);
  return ans;
}

//======================================================================

inline const string getClosingName(const string& key) {
  if ( ! isClosing(key) ) exit(0);
  int n = key.length();
  string ans = key.substr(2,n-3);
  return ans;
}

//#######################################################################

class Block {

private:

  string Name; // the name of the element
  string* Word; // the whole list of words to be processed
  string Open; // the key at which to start the process
  string Close; // the key at which to terminate the process
  int NB; // the number of subelements
  Array<Block> SubBlock;
  int NV; // the number of values
  Array<string> Value;

public:

  void initialize( string* word , const string& name = "" );
  void initialize( const string& FNAME , const string& name = "" );

  Block () {
    NB = 0; 
    NV = 0; 
  };

  Block (const string& FNAME, const string& BNAME="") {
    NB = 0; 
    NV = 0; 
    initialize(FNAME, BNAME);
  };

  ~Block () {
    //    printf("*** Destroying XML::Block (%s)\n", Name.c_str());
  };

  const int& NumberOfBlocks() const { return NB; };
  const int& NumberOfValues() const { return NV; };
  string getJoinedString();
  const string& getName() const;
  bool syntax_error();
  void read();
  int getInteger(const int i=0);
  double getDouble(const int i=0);
  string& getString(const int i=0);
  Block& getElement(const int i);
  Block& getElement(const string& name);
  Block& getElement(
    const string& name, const string& idkey, const int id );
  Block& operator[] (const int& i) { return SubBlock[i]; };
  Block& operator[] (const string& name);
  void dump(const string& prompt="");

};

//======================================================================

inline string Block::getJoinedString() {
  string s="";
  if ( NV == 0 ) return s;
  s = Value[0];
  for (int i=1; i<NV; i++) {
    s += " ";
    s += Value[i];
  }
  return s;
}

//======================================================================

inline const string& Block::getName() const { return Name; }

//======================================================================

inline void Block::initialize( string* word , const string& name ) {
//        printf("Block::initialize> Pass 1\n");
  Name = name;
  Word = word;
  if ( name == "" ) {
    Open = "";
    Close = "";
  } else {
    Open = "<"+name+">";
    Close = "</"+name+">";
  }
  read();
}

//======================================================================

void Block::initialize( const string& FNAME , const string& name ) {
//      printf("Block::initialize> Pass 2\n"); // koko

  FileReader F;
  F.open(FNAME.c_str());
  int NW;
  string* W;
  F.getWordList( NW , W ); //単語の数、最初の括弧<algorithm>
  initialize( W , name );
  delete [] W;
}

//======================================================================

inline bool Block::syntax_error() {
  NB = 0;
  NV = 0;
  int depth = 0;
  bool open = false;
  string SkipTo = "";
  int i=0;
  while(Word[i] != EOL) {
    string& w = Word[i++];
    if ( Name == "" ) {
      if ( isOpening(w) ) {
        const string name = getOpeningName(w);
        Name = name;
        Open = "<"+name+">";
        Close = "</"+name+">";
      }
    }
    if ( w == Open ) {
      depth++;
      open = true;
      continue;
    }
    if ( ! open ) continue;
    if ( w == Close ) {
      depth--;
      break;
    }
    if ( SkipTo != "" && w != SkipTo ) continue;
    if ( SkipTo != "" && w == SkipTo ) {
      depth--;
      SkipTo = "";
      continue;
    }
    if ( isCommentOpening(w) ) {
//cout << "  ... beginning of a comment" << endl;
      depth++;
      SkipTo = "-->";
      continue;
    }
    if ( isOpening(w) ) {
//cout << "  ... beginning of a subelement [" << w << "]" << endl;
      depth++;
      const string name = getOpeningName(w);
//cout << "opening name : " << name << endl;
      SkipTo = "</" + name + ">";
      NB++;
      continue;
    }
    if ( isClosing(w) ) {
      printf("Block::read> Error.\n");
      printf("  An unexpected closing tag %s\n", w.c_str());
      printf("  is detected in reading an element of name [%s].\n", 
        getName().c_str());
      return true;
    }
//cout << "  ... a value" << endl;
    NV++;
  }

  if ( depth != 0 ) {
    printf("Block::read> Error.\n");
    string expected = SkipTo;
    if ( expected == "" ) {
      expected = Close;
    }
    printf("  A missing closing tag %s\n", expected.c_str());
    printf("  is detected in reading an element of name [%s].\n", 
      getName().c_str());
    return true;
  }

  return false;

}

//======================================================================

void Block::read() {

//  printf("Block::read> Pass 1\n");
  if ( syntax_error() ) exit(0);

  if ( NV > 0 ) Value.init("Value",1,NV);
  if ( NB > 0 ) SubBlock.init("SubBlock",1,NB);

  bool open = false;
  string SkipTo = "";
  int ib = 0;
  int iv = 0;
  int i=0;
  while(true) {
    string& w = Word[i++];
//    cout << "### " << w << endl;
    if ( w == Open ) {
//printf("Opened. %s\n", w.c_str());
      open = true;
      continue;
    }
    if ( ! open ) continue;
    if ( w == Close ) break;
    if ( SkipTo != "" && w != SkipTo ) continue;
    if ( SkipTo != "" && w == SkipTo ) {
      SkipTo = "";
      continue;
    }
    if ( isCommentOpening(w) ) {
      SkipTo = "-->";
      continue;
    }
    if ( isOpening(w) ) {
      const string name = getOpeningName(w);
      SubBlock[ib].initialize(&w,name);
      ib++;
      SkipTo = "</" + name + ">";
      continue;
    }
    Value[iv] = w;
    iv++;
  }
}

//======================================================================

inline string& Block::getString(const int i) {
#ifdef DEB
  if ( i < 0 || i >= NV ) {
    printf("Block::getString> Error.\n");
    printf("  The argument (= %d) is out of the bounds [0,%d).\n", i, NV);
    exit(0);
  }
#endif
  return Value[i];
}

//======================================================================

inline int Block::getInteger(const int i) {
#ifdef DEB
  if ( i < 0 || i >= NV ) {
    printf("Block::getInteger> Error.\n");
    printf("  The argument (= %d) is out of the bounds [0,%d).\n", i, NV);
    exit(0);
  }
#endif
  string& s = getString(i);
  return atoi(s.c_str());
}

//======================================================================

inline double Block::getDouble(const int i) {
#ifdef DEB
  if ( i < 0 || i >= NV ) {
    printf("Block::getDouble> Error.\n");
    printf("  The argument (= %d) is out of the bounds [0,%d).\n", i, NV);
    exit(0);
  }
#endif
  string& s = getString(i);
  return (double)atof(s.c_str());
}

//======================================================================

inline Block& Block::getElement(const int i) {
  return SubBlock[i];
}

//======================================================================

inline Block& Block::getElement(const string& name) {
  int ib;
  for (ib=0; ib<NB; ib++) {
    if ( name == SubBlock[ib].getName() ) break;
  }
  if ( ib < NB ) { 
    return SubBlock[ib];
  } else {
    printf("Block::getElement> Error.\n");
    printf("  Tried to find a subblock '%s',\n",
      name.c_str());
    printf("  but no such block exists in a block '%s'.\n", 
      getName().c_str());
    exit(0);
  }
}

//======================================================================

inline Block& Block::getElement(
  const string& name, const string& idkey, const int id ) {
  int ib;
  for (ib=0; ib<NB; ib++) {
    if ( name == SubBlock[ib].getName() ) {
      Block& E = SubBlock[ib];
//    E.dump();
      int id0 = E[idkey].getInteger();
//    printf("id0 = %d, id = %d\n", id0, id);
      if ( id0 == id ) break;
    }
  }
  if ( ib < NB ) { 
    return SubBlock[ib];
  } else {
    printf("Block::getElement> Error.\n");
    printf("  Tried to find a subblock '%s',\n",
      name.c_str());
    printf("  but no such block exists in a block '%s'.\n", 
      getName().c_str());
    exit(0);
  }
}

//======================================================================

Block& Block::operator[] (const string& name) {
  Block& B = getElement(name);
  return B;
}

//======================================================================

inline void Block::dump(const string& prompt) {
  if ( NV == 1 && NB == 0 ) {
    printf("%s"     , prompt.c_str());
    printf("<%s>"   , Name.c_str());
    printf(" %s"    , Value[0].c_str());
    printf(" </%s>" , Name.c_str());
    printf("\n");
  } else {
    printf("%s<%s>\n", prompt.c_str(), Name.c_str());
    for (int i=0; i<NV; i++) {
      printf("%s  (%2d) %s\n", prompt.c_str(), i, Value[i].c_str());
    }
    for (int i=0; i<NB; i++) {
//    printf("%s  subelement %2d:\n", prompt.c_str(), i);
      SubBlock[i].dump(prompt+"  ");
    }
    printf("%s</%s>\n", prompt.c_str(), Name.c_str());
  }
}

//#######################################################################

} // 先頭の namespace XML { と対になる } ．

#endif
