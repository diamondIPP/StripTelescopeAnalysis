
#ifndef SaveToFile
#define SaveToFile

//
////void SaveCanvasPNG(TCanvas *canvas, char* location, char* file_name);
//void SaveCanvasPNG(TCanvas *canvas, string location, string file_name)
//{
//   char loc[500];
//   memcpy(loc,location.c_str(),strlen(location.c_str())+1);
//   char png[] = ".png";
//   char *pngloc = loc;
//
//   //Save .png file
//   strcat(pngloc,file_name.c_str());
//   strcat(pngloc,png);
//   char const *png_file = &pngloc[0]; //assigns address of first element of the file string to the char_file pointer
//   TImage *img = TImage::Create();
//   img->FromPad(canvas);
//   img->WriteImage(png_file);
//   cout << ".png file was created at: " << png_file << endl;
//   delete img;
//}
//
////void SaveCanvasC(TCanvas *canvas, char* location, char* file_name);
//void SaveCanvasC(TCanvas *canvas, string location, string file_name)
//{
//   char loc[500];
//   memcpy(loc,location.c_str(),strlen(location.c_str())+1);
//   char cmac[] = ".C";
//   char *cmacloc = loc;
//
//   //Saving .C macro
//   strcat(cmacloc,file_name.c_str());
//   strcat(cmacloc,cmac);
//   char const *cmac_file = &cmacloc[0];
//   canvas->SaveSource(cmac_file);
//   cout << ".c macro was created at: " << cmac_file << endl;
//}
//
////void SaveCanvasRoot(TCanvas *canvas, char* location, char* file_name);
//void SaveCanvasRoot(TCanvas *canvas, string location, string file_name)
//{
//    char loc[500];
//    memcpy(loc,location.c_str(),strlen(location.c_str())+1);
//    char rt[] = ".root";
//    char *rtloc = loc;
//    //Saving .root file
//    strcat(rtloc,file_name.c_str());
//    strcat(rtloc,rt);
//    char const *rt_file = &rtloc[0];
//    TObjArray list(0);
//    list.Add(canvas);
//    TFile f(rt_file,"recreate");
//    list.Write();
//    f.Close();
//    cout << ".root file was created at: " << rt_file << endl;
//}


#endif
