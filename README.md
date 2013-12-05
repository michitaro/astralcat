usage
-----
    > ls raw
    img1.cr2 img2.cr2
    
    # convert raw image into fits
    > mkdir fits
    > ./raw2fits -1 -o fits/img1.fits raw/img1.cr2
    > ./raw2fits -1 -o fits/img1.fits raw/img1.cr2
    
    # make catalog
    > mkdir catalog
    > ./sky --catalog=catalog/cat1.txt --detect='min_area=3 detect_threshold=4.5' --sky='type=localpoly' fits/img1.fits
    > ./sky --catalog=catalog/cat2.txt --detect='min_area=3 detect_threshold=4.5' --sky='type=localpoly' fits/img2.fits
    
    # stitch
    > ./stitch -o stack.fits -n5 catalog/* fits/*.fits



requirements
------------
  * C++0x compiler
    * http://gcc.gnu.org
      * CXX="g++ -std=gnu++11"
  * sfitsio
    * http://www.ir.isas.jaxa.jp/~cyamauch/sli/index.html
  * libraw
    * http://www.libraw.org
  * GSL
    * http://www.gnu.org/software/gsl/
  * boost
    * http://www.boost.org
