module mod_aero_coeff

      implicit none

CONTAINS

subroutine cpcrcm(APHIJ,EMIJ,  CLIFT, CDRAG_out, CMOME_out, ASLOP_out)

   real(kind=8)                        :: APHIJ, EMIJ
   real(kind=8), intent(out)           :: CLIFT
   real(kind=8), intent(out), optional :: CDRAG_out
   real(kind=8), intent(out), optional :: CMOME_out
   real(kind=8), intent(out), optional :: ASLOP_out

   real(kind=8) :: CDRAG
   real(kind=8) :: CMOME
   real(kind=8) :: ASLOP
   real(kind=8) :: SQT, C1, C2, C5, S, S2, S3, S4
   integer :: NEG
   real(kind=8), parameter :: pi = 3.1415926535



   NEG=1
   if (EMIJ>1.D0) EMIJ=.99D0
   SQT=sqrt(1.-EMIJ*EMIJ)
   !write(*,*) "SQT = ", SQT
   C1=1.-EMIJ
   C2=.22689*C1
   C5=EMIJ/(SQT*SQT)

   do while (APHIJ<0. .or. APHIJ>pi)

      if(APHIJ<0.) then
          APHIJ=-APHIJ
           NEG=-1*NEG
       endif
      if(APHIJ>pi) APHIJ=APHIJ-pi*2.
   end do

   if (APHIJ<C2) then
      ASLOP=5.7296/SQT
      CLIFT=ASLOP*APHIJ
      CDRAG=.006+.13131*APHIJ*APHIJ
      CMOME=1.4324*APHIJ/SQT
   else
      CDRAG=1.1233-.029894*cos(APHIJ)-1.00603*cos(2.*APHIJ)
      CDRAG=CDRAG+.003115*cos(3.*APHIJ)-.091487*cos(4.*APHIJ)
      CDRAG=CDRAG/SQT
      if (APHIJ<.34906) then
         CLIFT=.29269*C1+(1.3*EMIJ-.59)*APHIJ
         ASLOP = 1.3*EMIJ - 0.59
         C2=(.12217+.22689*EMIJ)*SQT
         CMOME=CLIFT/(4*C2)
         CLIFT=CLIFT/C2
      elseif (APHIJ<2.7402) then
         S=sin(APHIJ)
         S2=sin(2.*APHIJ)
         S3=sin(3.*APHIJ)
         S4=sin(4.*APHIJ)
         CLIFT=(.080373*S+1.04308*S2-.011059*S3+.023127*S4)/SQT
         ASLOP=(0.080373*cos(APHIJ) + 1.04308*2*cos(2.*APHIJ) - .011059*3*cos(3.*APHIJ) &
               + 0.023127*4*cos(4.*APHIJ))/SQT
         CMOME=(-.02827*S+.14022*S2-.00622*S3+.01012*S4)/SQT
      elseif (APHIJ<3.0020) then
         CLIFT=-(.4704+.10313*APHIJ)/SQT
         ASLOP = -0.10313/SQT
         CMOME=-(.4786+.02578*APHIJ)/SQT
      elseif(APHIJ<=pi)  then
         CLIFT=(-17.55+5.5864*APHIJ)/SQT
         ASLOP = 5.5864/SQT
         CMOME=(-12.5109+3.9824*APHIJ)/SQT
      endif
   endif

   if (NEG<=0) then
      CLIFT=-CLIFT
      CMOME=-CMOME
      APHIJ=-APHIJ
   endif
   CMOME=CMOME - .25*CLIFT

   if(present(CDRAG_out))  CDRAG_out = CDRAG
   if(present(CMOME_out))  CMOME_out = CMOME
   if(present(ASLOP_out))  ASLOP_out = ASLOP

end subroutine cpcrcm

end module mod_aero_coeff
