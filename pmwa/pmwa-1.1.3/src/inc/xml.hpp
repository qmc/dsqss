#ifndef XML_H
#define XML_H

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>


//io.h (here)##################
#define BLEN 256

inline int line_split(char* line, string* w) {
  string s(line);
  istringstream ist(s);
  int nw = 0;
  while ( ist >> w[nw++] );
  nw--;
  return nw;
}

class FileReader {

private:

  string EOL;
  char NAME[BLEN];
  char LINE[BLEN];
  string WORD[BLEN];
  int IL;
  int NW;
  ifstream INS;
  streampos top;

public:

  void open(const char* name) {
    strcpy(NAME,name);
    INS.open(NAME);
    top = INS.tellg();
    if (!INS) { printf("FILEIO::FILEIO> Error.\n"); exit(0); };
    IL = 0;
  };

  FileReader() {

    EOL="_E_O_L_";
  };

  void rewind() {
    INS.clear();
    INS.seekg(top);
  };

  int split() {
    NW = line_split( LINE, WORD);
    return NW;
  };

  bool read() {
    bool b = INS.getline( LINE, BLEN);
    if ( b ) { IL++; }
    return b;
  };

  string& word(int i) {
    if ( i<0 && i>=NW ) {
      printf("FileReader::word> Error.\n"); 
      exit(0);
    }
    return WORD[i]; 
  };


  void getWordList(int& NW, string*& W);

};


//io.h (end)################
//array.h (here)##################
//---------------------------------
class IndexSystem {

private:
  bool INI;
  string LBL; // label
  int  D;
  int* L;
  int  N;

public:
  void init( const int d, const int* l, const string& LBL0 = ""  );

  IndexSystem() { INI = false; };
  ~IndexSystem() { 
    if ( initialized() ) delete [] L;
  };

  bool initialized() const { return INI; };

};


//------------------------------
template<class C> class Array {

private:

  string LBL; // label
  int D;
  C* val;
  IndexSystem ID;
  void init(va_list& ap);

public:

  Array() { LBL = ""; val = 0; };
  ~Array();

  void reset();
  void init(const string& s, int d, ...);
  C& operator[] (const int i);

};


//array.h (end)##################


namespace XML {

  //class Block;


  class Block {

  private:

    string EOL;
    string Name; // the name of the element
    string* Word; // the whole list of words to be processed
    string Open; // the key at which to start the process
    string Close; // the key at which to terminate the process
    Array<Block> SubBlock;
    int NB; // the number of subelements
    int NV; // the number of values
    Array<string> Value;
    
  public:
    
    void initialize( string* word , const string& name = "" );
    void initialize( const string& FNAME , const string& name = "" );
    
    Block () {
      NB = 0; 
      NV = 0; 
      EOL="_E_O_L_";
    };

    ~Block () {
      //    printf("*** Destroying XML::Block (%s)\n", Name.c_str());
    };
    
    void read();
    ////
    bool syntax_error();
    ///
    
    Block& operator[] (const int& i);
    Block& operator[] (const string& name);
    
    const string& getName() const;
    int getInteger(const int i=0);
    double getDouble(const int i=0);
    string& getString(const int i=0);
    Block& getElement(const string& name);
    const int& NumberOfBlocks() const { return NB; };
    const int& NumberOfValues() const { return NV; };
    

  };


} // namespace XML { と対になる } ．


#endif
