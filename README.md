usage
-----
    ./raw2fits -o out.fits in.cr2
    ds9 &
    xpaset -p ds9 rgb
    xpaset -p ds9 rgb red   ; xpaset -p ds9 file out.fits'[0]'
    xpaset -p ds9 rgb green ; xpaset -p ds9 file out.fits'[1]'
    xpaset -p ds9 rgb blue  ; xpaset -p ds9 file out.fits'[2]'


requirements
------------
  * C++0x compiler
    * http://gcc.gnu.org
  * sfitsio
    * http://www.ir.isas.jaxa.jp/~cyamauch/sli/index.html
  * libraw
    * http://www.libraw.org
  * boost
    * http://www.boost.org
