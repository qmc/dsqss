
#ifndef IO_H
#define IO_H

//######################################################################

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "array.h"

#define BLEN 256

//######################################################################

using namespace std;

string EOL = "_E_O_L_";

//======================================================================

inline void reform_end_of_line(string& line) {
  int n = line.size();
  const char* a = line.c_str();
  if ( a[n-1] == 13 ) { // 13 stands for ^M
    line.replace(n-1,1,1,'\n');
  }
}

//======================================================================

inline int line_split(char* line, string* w) {
  string s(line);
  istringstream ist(s);
  int nw = 0;
  while ( ist >> w[nw++] );
  nw--;
  return nw;
}

//======================================================================

inline void get_line(FILE* F, char* line) {
  char buff[256];
  strcpy(buff,"");
  strcpy(line,"");
  fgets(buff,100,F);
  if (buff[0]==0 || buff[0]==13 || buff[0]==10) return;
  char* pch=strchr(buff,'#');
  if (pch != NULL) {
    int len=pch-buff;
    if (len > 0) {
      strncat(line,buff,len);
      strcat(line,"\n");
    }
  } else {
    strcat(line,buff);
  }
}

//======================================================================

inline void get_nbl(FILE* F, char* line) { // get the next non-blank line
  strcpy(line,"");
  while( ! strcmp(line,"") ) {
    if ( feof(F) ) break;
    get_line(F,line);
//  printf("%s\n",line);
  }
}

//======================================================================

inline int break_into_words(char* line, char* delim, char** word) { 
  char* last = line + strlen(line) - 1;
//printf( "line = '%s'\n", line);
//printf( "delimiter = '%s'\n", delim);
  char* w = line;
  int n = 0;
  while ( w != last ) {
    while ( w == strpbrk( w, delim) ) w++;
    char* p = strpbrk( w, delim);
    if ( p == 0 ) p = last;
//printf("\nw= %s", w);
//printf("p= %d, (*p) = '%c' (%d)\n", p, *p, *p);
    strncpy(word[n], w, p-w);
    strcat(word[n], "\0");
//printf("word[%d] = '%s'\n", n, word[n]);
    n++;
    w = p;
  }    
  return n;
}

//######################################################################

class FileReader {

private:

  char NAME[BLEN];
  char LINE[BLEN];
  int IL;
  int NW;
  string WORD[BLEN];
  ifstream INS;
  streampos top;
  streampos mark;

public:

  void open(const char* name) {
    strcpy(NAME,name);
    INS.open(NAME);
    top = INS.tellg();
    if (!INS) { printf("FileReader::open> Error. Unable to open %s\n", NAME); exit(0); };
    IL = 0;
  };

  FileReader() {};

  FileReader(char* name) {
    open(name);
  };

  void rewind() {
    INS.clear();
    INS.seekg(top);
  };

  void set_mark() {
    mark = INS.tellg();
  };

  void goto_mark() {
    INS.clear();
    INS.seekg(mark);
  };

  bool read() {
    bool b = INS.getline( LINE, BLEN);
    if ( b ) { IL++; }
    return b;
  };

  char* line() { return LINE; };

  int split() {
    NW = line_split( LINE, WORD);
    return NW;
  };

  string& word(int i) {
    if ( i<0 && i>=NW ) {
      printf("FileReader::word> Error.\n"); 
      exit(0);
    }
    return WORD[i]; 
  };

  int as_int(int i) {
    return atoi(WORD[i].c_str());
  };

  double as_double(int i) {
    return (double)atof(WORD[i].c_str());
  };

  void show() {
    cout << LINE << endl;
  };

  void dump() {
    cout << NAME << "[" << IL << "]> ";
    for (int i=0; i<NW; i++) {
      cout << " " << WORD[i];
    }
    cout << endl;
  };

  string get(char* key);

  int makeIndex( 
    const char* scope, const char* field, const char* k0, Array<int>& index );

  void getWordList(int& NW, string*& W);

};

//======================================================================

inline void FileReader::getWordList(int& NW, string*& W) {

  NW = 0;
  rewind();

  while(read()) NW += split(); 
  W = new string[NW+1];
  int iw = 0;
  rewind();
  while(read()) {
    int nw = split();
    for (int i=0; i<nw; i++) {
      W[iw] = word(i);
      iw++;
    }
  }
  W[NW] = EOL;

}

//======================================================================

inline string FileReader::get(char* key) {
  string ans;
  char key_open[BLEN];
  char key_close[BLEN];
  strcpy( key_open, "<");
  strcat( key_open, key);
  strcat( key_open, ">");
  strcpy( key_close, "</");
  strcat( key_close, key);
  strcat( key_close, ">");
  printf(" %s %s\n", key_open, key_close);
  rewind();
  bool active = false;
  while(read()) {
    char* l = line();
    //    int nw = split();
    string keyt = word(0);
    if (keyt == key_close) active = false;
    if (active) {
      string s = l;
      reform_end_of_line(s);
      ans += s;
    }
    if (keyt == key_open) active = true;
  }
  return ans;
}

//======================================================================

inline int FileReader::makeIndex( const char* scope, const char* field, const char* k0, Array<int>& index ) {

//  printf("FileReader::makeIndex> start.\n");

  string key(k0);
  string close = "</";
  close += scope;
  close += ">";
  string activate = "<";
  activate += field;
  activate += ">";
  string inactivate = "</";
  inactivate += field;
  inactivate += ">";

//  printf(" close= %s\n", close.c_str());
//  printf(" activate= %s\n", activate.c_str());
//  printf(" inactivate= %s\n", inactivate.c_str());

  int n=0;
  set_mark();
  bool active = false;
  while( read() ) {
//show();
    int nw = split();
    if ( nw == 0 ) continue;
    string k = word(0);
    if ( k == close ) break;
    if ( k == inactivate ) active = false;
    if ( active ) {
      if ( k == key ) {
//printf(" %s:%s:%s> %s\n", scope, field, k0, LINE);
        n++;
      }
    }
    if ( k == activate ) active = true;
  }

  index.init(1,n);

  int i=0;
  goto_mark();
  active = false;
  while(read()) {
    int nw = split();
    if ( nw == 0 ) continue;
    string k = word(0);
    if ( k == close ) break;
    if ( k == inactivate ) active = false;
    if (active) {
      if ( k == key ) {
        index[i] = as_int(1);
        i++;
      }
    }
    if ( k == activate ) active = true;
  }

  goto_mark();

//  printf("FileReader::makeIndex> end.\n");

  return n;
}

#endif
