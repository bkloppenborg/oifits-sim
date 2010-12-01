//============================================================================
// Name        : simulator.cpp
// Version     : 1.0
// Copyright ( C ) 2008-2009 Fabien Baron
//============================================================================

#include "simulator.h"
int main( )
{
  /*
   int nmag = 2 * 14;
   Row<double> values(nmag);
   for (int ii = 0; ii < nmag; ii++)
   values[ ii ] = phasemagrun((double) ii / 2.);
   cout << "Errors: " << values;
   ofstream file;
   file.open("phase_error_vs_mag.txt");
   file << values;
   file.close();
   */

  double test;
  test = phasemagrun(5.);
  return 0;

  /*
   double test;
   test= simplyfringes( 8. );
   return 0;
   */
}

double delaymagrun( double magnitude )
{

  // Observation in H band
  const double lowwav = 1.49e-6;
  const double hiwav = 1.78e-6;
  const double lambda0 = (lowwav + hiwav) / 2.;
  const double r0 = 0.5802;
  const double windspeed = 20.0; // windspeed in the dominant layer
  const double t0 = r0 / windspeed;

  const int nchannels = 4;
  const double frame_duration = t0 / 10.; // was t0/10 and worked fine

  // Basden
  const int n_coh = int(1.6 * t0 / frame_duration + .5);

  const int ndelta_coh = n_coh;

  const int n_incoh = int(30. * t0 / frame_duration + .5);

  const int nframes = n_coh * 50;
  cout << "Number of frames= " << nframes << " equivalent to " << nframes * frame_duration << " s" << endl;

  const double baseline = 50.0; //baseline in meters
  const double teldiameter = 1.4; // telescope diameter
  const int beamsize = 32;
  const int AO = 2;

  const double quantum_efficiency = 0.65;
  const double dark_current = 0.0;
  const double ADC_gain = 3.0;
  const double readnoise = 5.0;
  const int nreads = 1; // unused for the moment -- number of CCD hardware reads per integration, should reduce the RON

  // Integration
  const int ntimesteps = 1; // number of computed images accumulated per read (to simulate smearing if integration time is large)
  const int nwavesteps = 10;

  RanInit(random_number_seed, 1);

  // Setup target, stations and beams
  Source target = Source(NULL, 'H', magnitude, 2600., 14.7, 0.5, 0.0);

  ZernikeLib zernlib = ZernikeLib(3, beamsize);

  AtmosphereLayer layer1 = AtmosphereLayer(0, 1024, beamsize, lambda0, r0, (double) beamsize / teldiameter, windspeed, 0.0, 1e99, -1);
  AtmosphereLayer layer2 = AtmosphereLayer(0, 1024, beamsize, lambda0, r0, (double) beamsize / teldiameter, windspeed, 0.0, 1e99, -1);

  Station station1 = Station(0.0, 0.0, &layer1, 1.0, teldiameter);
  Station station2 = Station(baseline, 0.0, &layer2, 1.0, teldiameter);

  Array array1 = Array(0.0, 2, &station1, &station2 );

  Modulator modulator1 = Modulator(0, 0.0, 0.0, 0.0);
  Modulator modulator2 = Modulator(2, 2. * 6. * (double) n_coh * frame_duration, 60e-6, 0.0);

  double z0[ 3 ] = { 0.0, 0.0, 0.0 };
  //double z1[ 3 ] = { 0.0 , 0.0 , 0.0 };
  //double z2[ 3 ] = { 1.0 , 3e-2 , 0.0 };
  Delayline delayline1 = Delayline(0, z0, z0, 1e-3);
  Delayline delayline2 = Delayline(0, z0, z0, 1e-3);
  Spatialfilter spatialfilter1 = Spatialfilter(0, (double) beamsize / teldiameter, 10., 1e-3);

  SpectralMode wavelengths = SpectralMode(1, lowwav, hiwav, nchannels);

  Beam beam1 = Beam(&array1, 0,  &delayline1, &modulator1, &spatialfilter1, beamsize, AO, &zernlib);
  Beam beam2 = Beam(&array1, 1,  &delayline2, &modulator2, &spatialfilter1, beamsize, AO, &zernlib);

  Combiner combiner = Combiner(2.13e-6, 2, &beam1, &beam2);

  Detector CCD = Detector(nchannels, nreads, ntimesteps, nwavesteps, readnoise, quantum_efficiency, dark_current, ADC_gain);

  FTBasden FT = FTBasden(n_coh, ndelta_coh, n_incoh, &wavelengths, &modulator2, 1, 0);

  // i: frame index
  // t: time, t[i]  =  time for frame index i

  Row<double> t(nframes);
  Row<double> result(nframes / n_coh);
  Row<double> DLpos(nframes);
  Row<double> Atm(nframes);

  ofstream fringes, modulation;
  char outfile[ 60 ], modfile[ 60 ];
  sprintf(outfile, "fringes.txt");
  sprintf(modfile, "modulation.txt");
  fringes.open(outfile);
  modulation.open(modfile);

  for (int i = 0; i < nframes; i++)
  {
    t[ i ] = ((double) i + 1.) * frame_duration;
    CCD.FastOneReadP(t[ i ], &combiner, &wavelengths, &target);
    fringes << CCD.frame << "\n";
    modulation << beam2.modulator->GetDelay(t[ i ]) << "\n";
    delayline2.SetTargetPosition(0.);
    FT.StoreOneFrame(CCD.frame, t[ i ]);

    if (((i + 1) % n_coh) == 0)
    {
      cout << (i + 1) / n_coh - 1 << "\n";
      result[ (i + 1) / n_coh - 1 ] = FT.OPD_estimate;
      //delayline2.SetTargetPosition( - FT.OPD_estimate );
    }
    DLpos[ i ] = combiner.diagnostic_DL;
    Atm[ i ] = combiner.diagnostic_atm;
  }

  fringes.close();
  modulation.close();

  ofstream file1, file2;
  char str1[ 30 ], str2[ 30 ];
  sprintf(str1, "delay_DL_results%.1f.txt", (float) magnitude);
  sprintf(str2, "delay_atm_results%.1f.txt", (float) magnitude);
  file1.open(str1);
  file2.open(str2);
  file1 << DLpos;
  file2 << Atm;
  file1.close();
  file2.close();

  double mean_error = sqrt(variance(result));
  cout << "Magnitude: " << magnitude << " Delay Error: " << mean_error << endl;
  return mean_error;

}

double phasemagrun( double magnitude )
{
  // Initialize the RNG
  RanInit(random_number_seed, 1);

  // Spectrograph, observation in H band
  const double lowwav = 1.49e-6;
  const double hiwav = 1.78e-6;
  const double lambda0 = (lowwav + hiwav) / 2.;
  const int nchannels = 4;

  // Atmosphere
  const double r0 = 0.5802;
  const double windspeed = 10.0;
  const double t0 = r0 / windspeed;
  const double constant_piston = 0.0;

  // Telescope
  const double baseline = 50.0; //baseline in meters
  const double teldiameter = 1.4; // telescope diameter
  const int beamsize = 32;
  const int AO = 2;

  // Array
  const double quantum_efficiency = 0.65;
  const double dark_current = 0.0;
  const double ADC_gain = 1.0;
  const double readnoise = 5.0;
  const int nreads = 1; // unused for the moment -- number of CCD hardware reads per integration, should reduce the RON

  // Integration
  const int ntimesteps = 100; // number of computed images accumulated per read (to simulate smearing if integration time is large)
  const int nwavesteps = 10;
  // Cophasing
  const int nbins = 4; // used for phase tracking
  const int frames_per_bin = 1;
  const double frame_duration = (t0 / 80.) / double(frames_per_bin);
  const int nframes = 20 * nbins * frames_per_bin;

  // Setup target
  // type, band, magnitude, temperature, background mag/sq arcsecs, background aperture, declination, latitude
  Source target = Source(NULL, 'H', magnitude, 2600., 99., 0.5, 0.0);

  // Setup atmosphere, stations and beams
  ZernikeLib zernlib = ZernikeLib(AO + 1, beamsize);

  //	AtmosphereLayer::AtmosphereLayer(int generatorsize, int screensize, double lambda0, double r0,
  //	double pixellation, double windspeed, double winddirection, double outerscale, double scintillation_diameter)


  AtmosphereLayer layer1 = AtmosphereLayer(0, 1024, beamsize, lambda0, r0, (double) beamsize / teldiameter, windspeed, 0.0, 1e99, -1);
  AtmosphereLayer layer2 = AtmosphereLayer(0, 1024, beamsize, lambda0, r0, (double) beamsize / teldiameter, windspeed, 0.0, 1e99, -1);

  Station station1 = Station(0.0, 0.0, &layer1, 1.0, teldiameter);
  Station station2 = Station(baseline, 0.0, &layer2, 1.0, teldiameter);

  Array array1 = Array(0.0, 2, &station1, &station2 );

  SpectralMode wavelengths = SpectralMode(1, lowwav, hiwav, nchannels);
  cout << wavelengths.mean_wavenumber << endl;
  Modulator modulator1 = Modulator(0, 0.0, 0.0, 0.0);

  // Note: twice as many bins for ABCD over the period of a triangle, 4 for sawtooth
  Modulator modulator2 = Modulator(2, 2. * double(nbins * frames_per_bin) * frame_duration, 1. / wavelengths.mean_wavenumber[ 0 ], 0.0);

  double z1[ 3 ] =
  { 0.0, 0.0, 0.0 };
  double z2[ 3 ] =
  { 0.0, 0.0, 0.0 };
  Delayline delayline1 = Delayline(0, z1, z1, 1e-3);
  Delayline delayline2 = Delayline(0, z2, z2, 1e-3);
  Spatialfilter spatialfilter1 = Spatialfilter(0, (double) beamsize / teldiameter, 10., 1e-3);

  Beam beam1 = Beam(&array1, 0,  &delayline1, &modulator1, &spatialfilter1, beamsize, AO, &zernlib);
  Beam beam2 = Beam(&array1, 1,  &delayline2, &modulator2, &spatialfilter1, beamsize, AO, &zernlib);

  Combiner combiner = Combiner(constant_piston, 2, &beam1, &beam2);
  Detector CCD = Detector(nchannels, nreads, ntimesteps, nwavesteps, readnoise, quantum_efficiency, dark_current, ADC_gain);

  // Setup the FT algorithm
  FTphaseABCD2 FT = FTphaseABCD2(nbins, frames_per_bin, &wavelengths, &modulator2);
  //FTphaseABCD1 FT = FTphaseABCD1(nbins, frames_per_bin, &wavelengths, &modulator2);
  // FringeAcquisition Acq = FringeAcquisition( int type, )  // type 0 : coherencing 1: cophasing 2: data streaming
  // Setup Fringe acquisition module
  // - sets the timers/integration times
  // - get the frames whenever necessary, and passes them to the relevant FT module
  // - Postprocess phase ( servo ) before sending to delaylines


  // i: frame index
  // t: time, t[i]  =  time for frame index i

  //We get the result only after 4 bins
  Row<double> t(nframes);
  Row<double> result(nframes / (nbins * frames_per_bin));
  Row<double> DLpos(nframes);
  Row<double> Atm(nframes);
  ofstream fringes;
  char outfile[ 30 ];
  sprintf(outfile, "fringes.txt");
  fringes.open(outfile);

  cout << "\nMagnitude:" << magnitude << "\n";
  for (int i = 0; i < nframes; i++)
  {
    t[ i ] = ((double) i + 1.) * frame_duration;
    CCD.OneReadP(t[ i ], &combiner, &wavelengths, &target);
    fringes << CCD.frame;
    FT.StoreOneFrame(CCD.frame, t[ i ]);
    if (((i + 1) % (nbins * frames_per_bin)) == 0)
    {
      result[ (i + 1) / (nbins * frames_per_bin) - 1 ] = FT.opd_estimate;
      // delayline2.SetTargetPosition( - FT.phase_estimate );
      delayline2.SetTargetPosition(0);
    }
    DLpos[ i ] = combiner.diagnostic_DL;
    Atm[ i ] = combiner.diagnostic_atm;
  }

  fringes.close();

  double mean_error = sqrt(norm2(result));

  ofstream file1, file2;
  char str2[ 30 ]; //, str1[ 30 ];
  //sprintf( str1 , "OPD_DL_results%.1f.txt" , (float) magnitude );
  sprintf(str2, "OPD_atm_results%.1f.txt", (float) magnitude);
  //file1.open( str1 );
  file2.open(str2);
  //file1 << DLpos;
  file2 << Atm;
  //file1.close( );
  file2.close();
  return mean_error;

}

double simplyfringes( double magnitude )
{
  // Initialize the RNG
  RanInit(random_number_seed, 1);

  // Spectrograph, observation in H band
  const double lowwav = 1.49e-6;
  const double hiwav = 1.78e-6;
  const double lambda0 = (lowwav + hiwav) / 2.;
  const int nchannels = 1;

  // Atmosphere
  const double r0 = 0.5802;
  const double windspeed = 10.0;
  const double t0 = r0 / windspeed;
  const double constant_piston = 0.0;

  // Telescope
  const double teldiameter = 1.4; // telescope diameter
  const int beamsize = 32;
  const int AO = 2;

  // Array
  const double quantum_efficiency = 1.0;
  const double dark_current = 0.0;
  const double ADC_gain = 1.0;
  const double readnoise = 5.0;
  const int nreads = 1; // unused for the moment -- number of CCD hardware reads per integration, should reduce the RON

  // Integration
  const int ntimesteps = 1; // number of computed images accumulated per read (to simulate smearing if integration time is large)
  const int nwavesteps = 1;

  // Cophasing
  const double frame_duration = t0 / 20.;
  const int nframes = 3200;

  // Setup target
  // type, band, magnitude, temperature, background mag/sq arcsecs, background aperture, declination, latitude
  Source target = Source(NULL, 'H', magnitude, 2600., 99., 0.5, 0.0);

  // Setup atmosphere, stations and beams
  ZernikeLib zernlib = ZernikeLib(AO + 1, beamsize);

  AtmosphereLayer layer1 = AtmosphereLayer(0, 1024, beamsize, lambda0, r0, (double) beamsize / teldiameter, windspeed, 0.0, 1e99, -1);
  AtmosphereLayer layer2 = AtmosphereLayer(0, 1024, beamsize, lambda0, r0, (double) beamsize / teldiameter, windspeed, 0.0, 1e99, -1);
  AtmosphereLayer layer3 = AtmosphereLayer(0, 1024, beamsize, lambda0, r0, (double) beamsize / teldiameter, windspeed, 0.0, 1e99, -1);

  // Initialize the station positions to MROI ones                                                                                                          
  // West arm (indexes 0 to 9), North arm ( 10 to 18 ) and South arm ( 19 to 27 ) coordinates                                                               

  double MROI_north[ ] = 
    { 0.00, 1.89, 3.83, 5.66, 7.76, 11.73, 15.86, 23.96, 32.47, 49.05, 8.30, 10.82, 16.20, 21.37, 31.98, 42.38, 63.33, 83.99, 125.53,
     -7.16, -14.76, -21.97, -30.24, -45.83, -62.09, -93.93, -127.41, -192.64 };
  double MROI_east[ ] =
    { 0.00, -7.28, -14.94, -22.21, -30.55, -46.26, -62.64, -94.74, -128.48, -194.22, 2.07, 10.45, 15.67, 20.69, 31.00, 41.10, 61.44, 81.51,
      121.85, 2.04, 4.20, 6.24, 8.59, 13.01, 17.62, 26.66, 36.15, 54.66 };

  Station station1 = Station(MROI_east[ 0 ] , MROI_north[ 0 ], &layer1, 1.0, teldiameter);
  Station station2 = Station(MROI_east[ 10 ], MROI_north[ 10 ], &layer2, 1.0, teldiameter);
  Station station3 = Station(MROI_east[ 19 ], MROI_north[ 19 ], &layer3, 1.0, teldiameter);

  Array array1 = Array( 0.0 , 3, &station1, &station2, &station3 );

  SpectralMode wavelengths = SpectralMode(1, lowwav, hiwav, nchannels);
  cout << wavelengths.mean_wavenumber << endl;

  Modulator modulator1 = Modulator(4, 32. * frame_duration, 3. / wavelengths.mean_wavenumber[ 0 ], 0.0);
  Modulator modulator2 = Modulator(4, 32. * frame_duration, 2. / wavelengths.mean_wavenumber[ 0 ], 0.0);
  Modulator modulator3 = Modulator(4, 32. * frame_duration, -1. / wavelengths.mean_wavenumber[ 0 ], 0.0);

  double z1[ 3 ] =
  { 0.0, 0.0, 0.0 };
  double z2[ 3 ] =
  { 0.0, 0.0, 0.0 };
  double z3[ 3 ] =
  { 0.0, 0.0, 0.0 };

  Delayline delayline1 = Delayline(0, z1, z1, 1e-3);
  Delayline delayline2 = Delayline(0, z2, z2, 1e-3);
  Delayline delayline3 = Delayline(0, z3, z3, 1e-3);

  Spatialfilter spatialfilter1 = Spatialfilter(0, (double) beamsize / teldiameter, 10., 1e-3);

  Beam beam1 = Beam(&array1, 0,  &delayline1, &modulator1, &spatialfilter1, beamsize, AO, &zernlib);
  Beam beam2 = Beam(&array1, 1,  &delayline2, &modulator2, &spatialfilter1, beamsize, AO, &zernlib);
  Beam beam3 = Beam(&array1, 2,  &delayline3, &modulator3, &spatialfilter1, beamsize, AO, &zernlib);

  Combiner combiner = Combiner(constant_piston, 3, &beam1, &beam2, &beam3);
  Detector CCD = Detector(nchannels, nreads, ntimesteps, nwavesteps, readnoise, quantum_efficiency, dark_current, ADC_gain);

  // i: frame index
  // t: time, t[i]  =  time for frame index i

  //We get the result only after 4 bins
  Row<double> t(nframes);
  Row<double> DLpos(nframes);
  Row<double> Atm(nframes);
  ofstream fringes, modulation;
  char outfile[ 60 ], modfile[ 60 ];
  sprintf(outfile, "fringes.txt");
  sprintf(modfile, "modulation.txt");
  fringes.open(outfile);
  modulation.open(modfile);
  cout << "\nMagnitude:" << magnitude << "\n";
  for (int i = 0; i < nframes; i++)
  {
    t[ i ] = ((double) i + 1.) * frame_duration;
    CCD.OneReadP(t[ i ], &combiner, &wavelengths, &target);
    cout << "Frame " << i << "\r";
    fringes << CCD.frame << "\n";
    modulation << beam2.modulator->GetDelay(t[ i ]) << "\n";
    DLpos[ i ] = combiner.diagnostic_DL;
    Atm[ i ] = combiner.diagnostic_atm;
  }

  fringes.close();
  modulation.close();
  return 0.;
}

