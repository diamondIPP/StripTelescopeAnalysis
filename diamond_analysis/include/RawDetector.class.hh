/*
 * RawDetector.class.hh
 *
 *  Created on: 31.07.2011
 *      Author: Felix Bachmair
 */

#ifndef RAWDETECTOR_CLASS_HH_
#define RAWDETECTOR_CLASS_HH_

#include <iostream>
#include "TObject.h"
class RawDetector {
   public:
      //functions
      RawDetector();
      virtual ~RawDetector();
      void SetADC(int channel, unsigned short adc);
      unsigned short GetADC(Int_t channel);

   protected:
      unsigned short raw_detector_data[256];
};
#endif /* RAWDETECTOR_CLASS_HH_ */
