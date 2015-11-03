// the following fuctions are for graphical display
//                 (active with the compile option -D GRAPHIC)
//======================================================================

void Simulation::g_init() {
  int g_xorigin = 600;
  int g_yorigin = 10;
  int g_winwidth = 400;
  int g_winheight = 300;
  float g_winxmin = 0;
  float g_winxmax = 2.0;//(float)(L+1);
  float g_winymin = (float)(-0.1 * BETA);
  float g_winymax = (float)( 1.1 * BETA);
  if (D != 1) {
    printf(">>> Only 1D is supported for graphics.\n");
    exit(0);
  }
  if (L > 16) {
    printf(">>> L must not exceed 16 for graphics.\n");
    exit(0);
  }
  if (ALG.NX != 2) {
    printf(">>> Only 2-state model is supported for graphics.\n");
    exit(0);
  }
  printf(" %f %f %f %f\n",g_winxmin, g_winymin, g_winxmax, g_winymax);
  gsetinitialparsegeometry("%+d%+d",g_xorigin, g_yorigin);
  g_win = gopen(g_winwidth, g_winheight);
  winname(g_win, "test");
  window(g_win, g_winxmin, g_winymin, g_winxmax, g_winymax);
};  

//======================================================================

void Simulation::g_draw() {
  char* color[2] = {"blue", "red"};
  gclr(g_win);
  float x, x0, x1, y, y0, y1;
  for (int i=0; i<LAT.NSite; i++) { 
    ppS is = LAT.S[i]->Seg.begin();
    while(is != LAT.S[i]->Seg.end()) {
      x = (float)(LAT.S[i]->ID);
      y0 = (*((*is)->bottom))->tau;
      y1 = (*((*is)->top))->tau;
      if (y0 == BETA) y0 = 0.0;
      newcolor(g_win, color[(*is)->X]);
      line(g_win, x, y0, PSET);
      line(g_win, x, y1, PENDOWN);
      is++;
    }
  }
  newcolor(g_win, "white");
  for (int i=0; i<LAT.NBond; i++) {
    int is = LAT.B[i]->siteA;
    int js = LAT.B[i]->siteB;
    ppV iv = LAT.B[i]->Ver.begin();
    while(iv != LAT.B[i]->Ver.end()) {
      y = (*iv)->tau;
      if (is == 0 && js == L-1) {
        line(g_win,          0.5, y, PSET);
        line(g_win,          1.0, y, PENDOWN);
        line(g_win, (float)L    , y, PSET);
        line(g_win, (float)L+0.5, y, PENDOWN);
      } else {
        x0 = (float)(is+1);
        x1 = (float)(js+1);
        line(g_win, x0, y, PSET);
        line(g_win, x1, y, PENDOWN);
      }
      iv++;
    }
  }
  //  getchar();
};

//======================================================================
void Simulation::g_clear() {
  getchar();
  gclose(g_win);
};  
