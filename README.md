# pskons

This is a PSK31 program suitable for running in a Linux, MacOS, or
FreeBSD terminal (text) window.

It requires that fftw3, sndfile, and portaudio libraries already
be installed. To compile:

  make

To run, attaching to sound card number X:

  ./pskons -card X -out X

To find the numbers of sound card, run pskons without arguments.

You should see something like this in your terminal window:

`
  RX 1500  
  A 1605 0.9 9xc     How copy? BTU Fred Smith, 3D6P de CE1AQ pse kn   
  B  812 0.6 ear Samantha t r efort (RSQ) : 5 599 EFD*:fB  S/uomy dB1N   b  
  C 1618 0.7 elieTelo  
  D 1621 0.7 elieTelo  
  E 1165 0.6 oai copy  GD DX Luck and Health  EQSL only TUqs n qsb 73     !e  
  F 1181 0.6     atC5te ioTee  ee eeaea e t etleeea e  PF i  l iei  a eMoa   
  G 1590 0.8 i i e e ieeleetTeei%lte %lae e i-eate  e%a aVa  Pe elieTele  
  H 1566 0.6    t re   r      r  H  ct ttr r    t r D    t   ret i   ne ,ae i  
  I 1637 0.6  a  ta   T      e  eee   e   i     t      e        ta  rt  t, e   
  J 1509 0.5   
  K 2129 0.6  eeo e eee-  T soel eo1eute   rt  ar    het  
  -   
  -   
  -   
  -   
  -   
  -   
  >   
  >   
  >   
  >   
  >   
`

The lines starting with A-K are a Digipan-like display of the most
plausible signals (typically mostly not signals at all). You can
select one of these signals to talk to with control-A followed by the
(lower case) letter of the signal. Then the lines starting with "-"
will show received text from that signal. Send by typing; your text
will show up on the lines starting with ">". Type control-C or
control-X to quit. The file qso-trace.txt will contain the text of
signals you have selected with control-A.
