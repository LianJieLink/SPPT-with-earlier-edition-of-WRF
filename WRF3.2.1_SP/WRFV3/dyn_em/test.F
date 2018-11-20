subroutine UPDATE_STOCH(                                         & 
                       SPSTREAMFORCS,SPSTREAMFORCC,SPTFORCS,SPTFORCC,  &
                       SPT_AMP,SPSTREAM_AMP,                           &
                       ids, ide, jds, jde, kds, kde,                   &
                       ims, ime, jms, jme, kms, kme,                   &
                       its, ite, jts, jte, kts, kte                    )

      IMPLICIT NONE

      REAL, DIMENSION( ids:ide,jds:jde)      :: ZRANDNOSS1,ZRANDNOSC1,ZRANDNOSS2,ZRANDNOSC2
      REAL, DIMENSION (ims:ime,jms:jme)      :: SPSTREAMFORCS,SPSTREAMFORCC,SPTFORCS,SPTFORCc,SPSTREAM_AMP,SPT_AMP
      INTEGER , INTENT(IN)     ::               ids, ide, jds, jde, kds, kde,   &
                                                ims, ime, jms, jme, kms, kme,   &
                                                its, ite, jts, jte, kts, kte
    
      REAL :: Z
      REAL, PARAMETER :: thresh = 3.0 
      INTEGER ::IL, IK,LMAX,KMAX
      LOGICAL :: LGAUSS

      KMAX=(jde-jds)+1 !NLAT
      LMAX=(ide-ids)+1 !NATX
!     Pick the distribution of the noise
!     Random noise uses global indexes to ensure necessary symmetries and anti-symmetries
!     of random forcing when run on multiple processors
      LGAUSS=.true.
      IF (LGAUSS) then
        DO IK=jds,jde
          DO IL=ids,ide
           do
            call gauss_noise(z)
            if (abs(z)<thresh) exit
           ENDDO
           ZRANDNOSS1(IL,IK)=z
           do
            call gauss_noise(z)
            if (abs(z)<thresh) exit
           ENDDO
           ZRANDNOSC1(IL,IK)=z
           do
            call gauss_noise(z)
            if (abs(z)<thresh) exit
           ENDDO
           ZRANDNOSS2(IL,IK)=z
           do
            call gauss_noise(z)
            if (abs(z)<thresh) exit
           ENDDO
           ZRANDNOSC2(IL,IK)=z
          ENDDO
        ENDDO
      ELSE
        DO IK=jds,jde
          DO IL=ids,ide
            CALL RANDOM_NUMBER(z)
            ZRANDNOSS1(IL,IK)=z-0.5
            CALL RANDOM_NUMBER(z)
            ZRANDNOSC1(IL,IK)=z-0.5
            CALL RANDOM_NUMBER(z)
            ZRANDNOSS2(IL,IK)=z-0.5
            CALL RANDOM_NUMBER(z)
            ZRANDNOSC2(IL,IK)=z-0.5
          ENDDO
        ENDDO
      ENDIF

!     Note: There are symmetries and anti-symmetries to ensure real-valued back transforms
! for symmetric part: left and right half axis symmetric

      DO IK=jts,jte
      if ((IK.le.(KMAX/2+1)) .and. (IK>1)) then ! Upper half
        DO IL=its,ite
          SPSTREAMFORCC(IL,IK)  = (1.-ALPH)*SPSTREAMFORCC(IL,IK) + SPSTREAM_AMP(IL,IK)* ZRANDNOSC1(IL,IK)
          SPSTREAMFORCS(IL,IK)  = (1.-ALPH)*SPSTREAMFORCS(IL,IK) + SPSTREAM_AMP(IL,IK)* ZRANDNOSS1(IL,IK) 
          SPTFORCC(IL,IK)       = (1.-ALPH)*SPTFORCC(IL,IK)      + SPT_AMP(IL,IK)     * ZRANDNOSC2(IL,IK)  
          SPTFORCS(IL,IK)       = (1.-ALPH)*SPTFORCS(IL,IK)      + SPT_AMP(IL,IK)     * ZRANDNOSS2(IL,IK)  
        ENDDO
      ELSEIF (IK==1) then
        DO IL=its,ite
        if ((IL.le.(LMAX/2+1))) then
          SPSTREAMFORCC(IL,IK)  = (1.-ALPH)*SPSTREAMFORCC(IL,IK) + SPSTREAM_AMP(IL,IK)* ZRANDNOSC1(IL,IK)
          SPSTREAMFORCS(IL,IK)  = (1.-ALPH)*SPSTREAMFORCS(IL,IK) + SPSTREAM_AMP(IL,IK)* ZRANDNOSS1(IL,IK) 
          SPTFORCC(IL,IK)       = (1.-ALPH)*SPTFORCC(IL,IK)      + SPT_AMP(IL,IK)     * ZRANDNOSC2(IL,IK)  
          SPTFORCS(IL,IK)       = (1.-ALPH)*SPTFORCS(IL,IK)      + SPT_AMP(IL,IK)     * ZRANDNOSS2(IL,IK)  
        elseif ((IL.gt.(LMAX/2+1))) then
          SPSTREAMFORCC(IL,IK)  = (1.-ALPH)*SPSTREAMFORCC(IL,IK) + SPSTREAM_AMP(IL,IK)* ZRANDNOSC1(LMAX-IL+2,IK)
          SPSTREAMFORCS(IL,IK)  = (1.-ALPH)*SPSTREAMFORCS(IL,IK) - SPSTREAM_AMP(IL,IK)* ZRANDNOSS1(LMAX-IL+2,IK) 
          SPTFORCC(IL,IK)       = (1.-ALPH)*SPTFORCC(IL,IK)      + SPT_AMP(IL,IK)     * ZRANDNOSC2(LMAX-IL+2,IK)  
          SPTFORCS(IL,IK)       = (1.-ALPH)*SPTFORCS(IL,IK)      - SPT_AMP(IL,IK)     * ZRANDNOSS2(LMAX-IL+2,IK)  
        endif
        ENDDO
      ENDIF
      ENDDO

      DO IK=jts,jte
      if (IK.gt.(KMAX/2+1)) then ! Lower half
        DO IL=its,ite
          if (IL.le.(LMAX/2+1).and.(IL.gt.1)) then !lower left 
           SPSTREAMFORCC(IL,IK) = (1.-ALPH)* SPSTREAMFORCC(IL,IK) + SPSTREAM_AMP(IL,IK) * ZRANDNOSC1(LMAX-IL+2,KMAX-IK+2)
           SPSTREAMFORCS(IL,IK) = (1.-ALPH)* SPSTREAMFORCS(IL,IK) - SPSTREAM_AMP(IL,IK) * ZRANDNOSS1(LMAX-IL+2,KMAX-IK+2)
           SPTFORCC(IL,IK)      = (1.-ALPH)* SPTFORCC(IL,IK)      + SPT_AMP(IL,IK)      * ZRANDNOSC2(LMAX-IL+2,KMAX-IK+2)
           SPTFORCS(IL,IK)      = (1.-ALPH)* SPTFORCS(IL,IK)      - SPT_AMP(IL,IK)      * ZRANDNOSS2(LMAX-IL+2,KMAX-IK+2)
          elseif (IL.eq.1) then !don't exceed index
           SPSTREAMFORCC(IL,IK) = (1.-ALPH)* SPSTREAMFORCC(IL,IK) + SPSTREAM_AMP(IL,IK) * ZRANDNOSC1(        1,KMAX-IK+2)
           SPSTREAMFORCS(IL,IK) = (1.-ALPH)* SPSTREAMFORCS(IL,IK) - SPSTREAM_AMP(IL,IK) * ZRANDNOSS1(        1,KMAX-IK+2)
           SPTFORCC(IL,IK)      = (1.-ALPH)* SPTFORCC(IL,IK)      + SPT_AMP(IL,IK)      * ZRANDNOSC2(        1,KMAX-IK+2)
           SPTFORCS(IL,IK)      = (1.-ALPH)* SPTFORCS(IL,IK)      - SPT_AMP(IL,IK)      * ZRANDNOSS2(        1,KMAX-IK+2)
          elseif (IL.gt.(LMAX/2+1)) then !lower right
           SPSTREAMFORCC(IL,IK) = (1.-ALPH)* SPSTREAMFORCC(IL,IK) + SPSTREAM_AMP(IL,IK) * ZRANDNOSC1(LMAX-IL+2,KMAX-IK+2)
           SPSTREAMFORCS(IL,IK) = (1.-ALPH)* SPSTREAMFORCS(IL,IK) - SPSTREAM_AMP(IL,IK) * ZRANDNOSS1(LMAX-IL+2,KMAX-IK+2)
           SPTFORCC(IL,IK)      = (1.-ALPH)* SPTFORCC(IL,IK)      + SPT_AMP(IL,IK)      * ZRANDNOSC2(LMAX-IL+2,KMAX-IK+2)
           SPTFORCS(IL,IK)      = (1.-ALPH)* SPTFORCS(IL,IK)      - SPT_AMP(IL,IK)      * ZRANDNOSS2(LMAX-IL+2,KMAX-IK+2)
          endif
        ENDDO
      endif
      ENDDO


     END subroutine UPDATE_STOCH
