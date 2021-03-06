!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine prints out the current value of variables at all specified
!   time series locations that are within the current patch.
!
! Michael G. Duda -- 25 August 2005
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calc_ts_locations( grid )

   USE module_domain, ONLY : domain, get_ijk_from_grid
   USE module_configure, ONLY : model_config_rec, grid_config_rec_type, model_to_grid_config_rec
   USE module_dm, ONLY : wrf_dm_min_real
   USE module_llxy
   USE module_state_description

   IMPLICIT NONE

   ! Arguments
   TYPE (domain), INTENT(INOUT) :: grid

   ! Externals
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
   INTEGER, EXTERNAL :: get_unused_unit

   ! Local variables
   INTEGER :: ntsloc_temp
   INTEGER :: i, k, iunit
   REAL :: ts_rx, ts_ry, ts_xlat, ts_xlong, ts_hgt
   REAL :: known_lat, known_lon
   CHARACTER (LEN=132) :: message
   TYPE (PROJ_INFO) :: ts_proj
   TYPE (grid_config_rec_type) :: config_flags

   INTEGER :: ids, ide, jds, jde, kds, kde,        &
              ims, ime, jms, jme, kms, kme,        &
              ips, ipe, jps, jpe, kps, kpe,        &
              imsx, imex, jmsx, jmex, kmsx, kmex,  &
              ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
              imsy, imey, jmsy, jmey, kmsy, kmey,  &
              ipsy, ipey, jpsy, jpey, kpsy, kpey


   IF ( grid%ntsloc .LE. 0 ) RETURN

#if ((EM_CORE == 1) && (DA_CORE != 1))
   IF ( grid%dfi_stage == DFI_FST ) THEN
#endif
      CALL get_ijk_from_grid ( grid ,                               &
                               ids, ide, jds, jde, kds, kde,        &
                               ims, ime, jms, jme, kms, kme,        &
                               ips, ipe, jps, jpe, kps, kpe,        &
                               imsx, imex, jmsx, jmex, kmsx, kmex,  &
                               ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
                               imsy, imey, jmsy, jmey, kmsy, kmey,  &
                               ipsy, ipey, jpsy, jpey, kpsy, kpey )
   
      CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )
   
      ! Set up map transformation structure
      CALL map_init(ts_proj)
   
      IF (ips <= 1 .AND. 1 <= ipe .AND. &
          jps <= 1 .AND. 1 <= jpe) THEN
         known_lat = grid%xlat(1,1)
         known_lon = grid%xlong(1,1)
      ELSE
         known_lat = 9999.
         known_lon = 9999.
      END IF
      known_lat = wrf_dm_min_real(known_lat)
      known_lon = wrf_dm_min_real(known_lon)
   
      ! Mercator
      IF (config_flags%map_proj == PROJ_MERC) THEN
         CALL map_set(PROJ_MERC, ts_proj,               &
                      truelat1 = config_flags%truelat1, &
                      lat1     = known_lat,             &
                      lon1     = known_lon,             &
                      knowni   = 1.,                    &
                      knownj   = 1.,                    &
                      dx       = config_flags%dx)
   
      ! Lambert conformal
      ELSE IF (config_flags%map_proj == PROJ_LC) THEN
         CALL map_set(PROJ_LC, ts_proj,                  &
                      truelat1 = config_flags%truelat1,  &
                      truelat2 = config_flags%truelat2,  &
                      stdlon   = config_flags%stand_lon, &
                      lat1     = known_lat,              &
                      lon1     = known_lon,              &
                      knowni   = 1.,                     &
                      knownj   = 1.,                     &
                      dx       = config_flags%dx)
   
      ! Polar stereographic
      ELSE IF (config_flags%map_proj == PROJ_PS) THEN
         CALL map_set(PROJ_PS, ts_proj,                  &
                      truelat1 = config_flags%truelat1,  &
                      stdlon   = config_flags%stand_lon, &
                      lat1     = known_lat,              &
                      lon1     = known_lon,              &
                      knowni   = 1.,                     &
                      knownj   = 1.,                     &
                      dx       = config_flags%dx)
   
#if (EM_CORE == 1)
      ! Cassini (global ARW)
      ELSE IF (config_flags%map_proj == PROJ_CASSINI) THEN
         CALL map_set(PROJ_CASSINI, ts_proj,                            &
                      latinc   = grid%dy*360.0/(2.0*EARTH_RADIUS_M*PI), &
                      loninc   = grid%dx*360.0/(2.0*EARTH_RADIUS_M*PI), & 
                      lat1     = known_lat,                             &
                      lon1     = known_lon,                             &
                      lat0     = config_flags%pole_lat,                 &
                      lon0     = config_flags%pole_lon,                 &
                      knowni   = 1.,                                    &
                      knownj   = 1.,                                    &
                      stdlon   = config_flags%stand_lon)
#endif

      ! Rotated latitude-longitude
      ELSE IF (config_flags%map_proj == PROJ_ROTLL) THEN
         CALL map_set(PROJ_ROTLL, ts_proj,                      &
! I have no idea how this should work for NMM nested domains
                      ixdim    = grid%e_we-1,                   &
                      jydim    = grid%e_sn-1,                   &
                      phi      = real(grid%e_sn-2)*grid%dy/2.0, &
                      lambda   = real(grid%e_we-2)*grid%dx,     &
                      lat1     = config_flags%cen_lat,          &
                      lon1     = config_flags%cen_lon,          &
                      latinc   = grid%dy,                       &
                      loninc   = grid%dx,                       &
                      stagger  = HH)
   
      END IF
   
      ! Determine time series locations for domain
      IF (.NOT. grid%have_calculated_tslocs) THEN
         grid%have_calculated_tslocs = .TRUE.
         WRITE(message, '(A43,I3)') 'Computing time series locations for domain ', grid%id
         CALL wrf_message(message)
   
         ntsloc_temp = 0
         DO k=1,grid%ntsloc
   
            IF (config_flags%map_proj == 0) THEN  ! For idealized cases, no map transformation needed
               ts_rx = grid%lattsloc(k)           ! NB: (x,y) = (lat,lon) rather than (x,y) = (lon,lat)
               ts_ry = grid%lontsloc(k)
            ELSE
               CALL latlon_to_ij(ts_proj, grid%lattsloc(k), grid%lontsloc(k), ts_rx, ts_ry)
            END IF
            

            ntsloc_temp = ntsloc_temp + 1
            grid%itsloc(ntsloc_temp) = NINT(ts_rx)
            grid%jtsloc(ntsloc_temp) = NINT(ts_ry)
            grid%id_tsloc(ntsloc_temp) = k
   
            ! Is point outside of domain (or on the edge of domain)?
            IF (grid%itsloc(ntsloc_temp) < ids .OR. grid%itsloc(ntsloc_temp) > ide .OR. &
                grid%jtsloc(ntsloc_temp) < jds .OR. grid%jtsloc(ntsloc_temp) > jde) THEN
               ntsloc_temp = ntsloc_temp - 1
   
            END IF
   
         END DO
   
         grid%next_ts_time = 1
   
         grid%ntsloc_domain = ntsloc_temp
   
         DO k=1,grid%ntsloc_domain
   
            ! If location is outside of patch, we need to get lat/lon of TS grid cell from another patch
            IF (grid%itsloc(k) < ips .OR. grid%itsloc(k) > ipe .OR. &
                grid%jtsloc(k) < jps .OR. grid%jtsloc(k) > jpe) THEN
               ts_xlat  = 1.E30
               ts_xlong = 1.E30
               ts_hgt   = 1.E30
            ELSE
               ts_xlat  = grid%xlat(grid%itsloc(k),grid%jtsloc(k))
               ts_xlong = grid%xlong(grid%itsloc(k),grid%jtsloc(k))
#if (EM_CORE == 1)
               ts_hgt   = grid%ht(grid%itsloc(k),grid%jtsloc(k))
#endif
            END IF
#if DM_PARALLEL
            ts_xlat  = wrf_dm_min_real(ts_xlat)
            ts_xlong = wrf_dm_min_real(ts_xlong)
            ts_hgt   = wrf_dm_min_real(ts_hgt)
#endif
   
            IF ( wrf_dm_on_monitor() ) THEN

               iunit = get_unused_unit()
               IF ( iunit <= 0 ) THEN
                  CALL wrf_error_fatal('Error in calc_ts_locations: could not find a free Fortran unit.')
               END IF

               WRITE(grid%ts_filename(k),'(A)') TRIM(grid%nametsloc(grid%id_tsloc(k)))//'.d00.TS'
               i = LEN_TRIM(grid%ts_filename(k))
               WRITE(grid%ts_filename(k)(i-4:i-3),'(I2.2)') grid%id
               OPEN(UNIT=iunit, FILE=TRIM(grid%ts_filename(k)), FORM='FORMATTED', STATUS='REPLACE')
#if (EM_CORE == 1)
               WRITE(UNIT=iunit, &
                     FMT='(A26,I2,I3,A6,A2,F7.3,A1,F8.3,A3,I4,A1,I4,A3,F7.3,A1,F8.3,A2,F6.1,A7)') &
                     grid%desctsloc(grid%id_tsloc(k))//' ', grid%id, grid%id_tsloc(k), &
                     ' '//grid%nametsloc(grid%id_tsloc(k)), &
                     ' (', grid%lattsloc(grid%id_tsloc(k)), ',', grid%lontsloc(grid%id_tsloc(k)), ') (', &
                     grid%itsloc(k), ',', grid%jtsloc(k), ') (', &
                     ts_xlat, ',', ts_xlong, ') ', &
                     ts_hgt,' meters'
#else
               WRITE(UNIT=iunit, &
                     FMT='(A26,I2,I3,A6,A2,F7.3,A1,F8.3,A3,I4,A1,I4,A3,F7.3,A1,F8.3,A2)') &
                     grid%desctsloc(grid%id_tsloc(k))//' ', grid%id, grid%id_tsloc(k), &
                     ' '//grid%nametsloc(grid%id_tsloc(k)), &
                     ' (', grid%lattsloc(grid%id_tsloc(k)), ',', grid%lontsloc(grid%id_tsloc(k)), ') (', &
                     grid%itsloc(k), ',', grid%jtsloc(k), ') (', &
                     ts_xlat, ',', ts_xlong, ') '
#endif
               CLOSE(UNIT=iunit)
            END IF
         END DO
   
      END IF
#if ((EM_CORE == 1) && (DA_CORE != 1))
   END IF
#endif

END SUBROUTINE calc_ts_locations


SUBROUTINE calc_ts( grid )

   USE module_domain
   USE module_model_constants

   IMPLICIT NONE

   ! Arguments
   TYPE (domain), INTENT(INOUT) :: grid

   LOGICAL, EXTERNAL :: wrf_dm_on_monitor

   ! Local variables
   INTEGER :: i, k, mm, n, ix, iy, rc
   REAL :: earth_u, earth_v, output_t, output_q, clw, xtime_minutes
   REAL, ALLOCATABLE, DIMENSION(:) :: p8w

   ! Parameter ts_model_level:  
       ! TRUE to output T, Q, and wind at lowest model level
       ! FALSE to output T and Q at 2-m and wind at 10-m diagnostic levels:
   LOGICAL, PARAMETER :: ts_model_level = .FALSE.  

   IF ( grid%ntsloc_domain .LE. 0 ) RETURN

#if ((EM_CORE == 1) && (DA_CORE != 1))
   IF ( grid%dfi_opt /= DFI_NODFI .AND. grid%dfi_stage /= DFI_FST ) RETURN
#endif

   n = grid%next_ts_time

   ALLOCATE(p8w(grid%sm32:grid%em32))

   DO i=1,grid%ntsloc_domain

      ix = grid%itsloc(i)
      iy = grid%jtsloc(i)
  
      IF (grid%sp31 <= ix .AND. ix <= grid%ep31 .AND. &
          grid%sp33 <= iy .AND. iy <= grid%ep33) THEN
       
         IF (ts_model_level) THEN
   
            !
            ! Output from the lowest model computational level:
            !
#if (EM_CORE == 1)
            earth_u = grid%u_2(ix,1,iy)*grid%cosa(ix,iy)-grid%v_2(ix,1,iy)*grid%sina(ix,iy)
            earth_v = grid%v_2(ix,1,iy)*grid%cosa(ix,iy)+grid%u_2(ix,1,iy)*grid%sina(ix,iy)
            output_t = grid%t_2(ix,1,iy)
#else
            earth_u = grid%u(ix,1,iy)
            earth_v = grid%v(ix,1,iy)
            output_t = grid%t(ix,1,iy)
#endif
            output_q = grid%moist(ix,1,iy,P_QV)
   
         ELSE
   
            !
            ! Output at 2-m and 10-m diagnostic levels:
            !
#if (EM_CORE == 1)
            earth_u = grid%u10(ix,iy)*grid%cosa(ix,iy)-grid%v10(ix,iy)*grid%sina(ix,iy)
            earth_v = grid%v10(ix,iy)*grid%cosa(ix,iy)+grid%u10(ix,iy)*grid%sina(ix,iy)
            output_q = grid%q2(ix,iy)
#else
            earth_u = grid%u10(ix,iy)
            earth_v = grid%v10(ix,iy)
            output_q = grid%qsfc(ix,iy)
#endif
            output_t = grid%t2(ix,iy)
   
         END IF
   
#if (EM_CORE == 1)
         ! Calculate column-integrated liquid/ice  (kg/m^2 or mm)
         CALL calc_p8w(grid, ix, iy, p8w, grid%sm32, grid%em32)
         clw=0.
         DO mm = 1, num_moist
            IF ( (mm == P_QC) .OR. (mm == P_QR) .OR. (mm == P_QI) .OR. &
                 (mm == P_QS) .OR. (mm == P_QG) ) THEN
               DO k=grid%sm32,grid%em32-1
                  clw=clw+grid%moist(ix,k,iy,mm)*(p8w(k)-p8w(k+1))
               END DO
           END IF
         END DO
         clw = clw / g
#endif
   
         CALL domain_clock_get( grid, minutesSinceSimulationStart=xtime_minutes )
         grid%ts_hour(n,i) = xtime_minutes / 60.
         grid%ts_u(n,i)    = earth_u
         grid%ts_v(n,i)    = earth_v
         grid%ts_t(n,i)    = output_t
         grid%ts_q(n,i)    = output_q
         grid%ts_psfc(n,i) = grid%psfc(ix,iy)
#if (EM_CORE == 1)
         grid%ts_glw(n,i)  = grid%glw(ix,iy)
         grid%ts_gsw(n,i)  = grid%gsw(ix,iy)
         grid%ts_hfx(n,i)  = grid%hfx(ix,iy)
         grid%ts_lh(n,i)   = grid%lh(ix,iy)
         grid%ts_clw(n,i)  = clw
         grid%ts_rainc(n,i)  = grid%rainc(ix,iy)
         grid%ts_rainnc(n,i) = grid%rainnc(ix,iy)
         grid%ts_tsk(n,i)  = grid%tsk(ix,iy)
#else
         grid%ts_tsk(n,i)  = grid%nmm_tsk(ix,iy)
#endif
         grid%ts_tslb(n,i) = grid%tslb(ix,1,iy)
   
      ELSE
  
         grid%ts_hour(n,i) = 1.E30
         grid%ts_u(n,i)    = 1.E30
         grid%ts_v(n,i)    = 1.E30
         grid%ts_t(n,i)    = 1.E30
         grid%ts_q(n,i)    = 1.E30
         grid%ts_psfc(n,i) = 1.E30
#if (EM_CORE == 1)
         grid%ts_glw(n,i)  = 1.E30
         grid%ts_gsw(n,i)  = 1.E30
         grid%ts_hfx(n,i)  = 1.E30
         grid%ts_lh(n,i)   = 1.E30
         grid%ts_clw(n,i)  = 1.E30
         grid%ts_rainc(n,i)  = 1.E30
         grid%ts_rainnc(n,i) = 1.E30
#endif
         grid%ts_tsk(n,i)  = 1.E30
         grid%ts_tslb(n,i) = 1.E30
   
      END IF
   END DO

   DEALLOCATE(p8w)
 
   grid%next_ts_time = grid%next_ts_time + 1

   IF ( grid%next_ts_time > grid%ts_buf_size ) CALL write_ts(grid)

END SUBROUTINE calc_ts


SUBROUTINE write_ts( grid )

   USE module_domain, ONLY : domain
   USE module_dm, ONLY : wrf_dm_min_reals
   USE module_state_description

   IMPLICIT NONE

   ! Arguments
   TYPE (domain), INTENT(INOUT) :: grid

   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
   INTEGER, EXTERNAL :: get_unused_unit

   ! Local variables
   INTEGER :: i, n, ix, iy, iunit
   REAL, ALLOCATABLE, DIMENSION(:,:) :: ts_buf

   IF ( grid%ntsloc_domain .LE. 0 ) RETURN

#if ((EM_CORE == 1) && (DA_CORE != 1))
   IF ( grid%dfi_opt /= DFI_NODFI .AND. grid%dfi_stage /= DFI_FST ) RETURN
#endif

#ifdef DM_PARALLEL
   ALLOCATE(ts_buf(grid%ts_buf_size,grid%max_ts_locs))

   ts_buf(:,:) = grid%ts_hour(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_hour(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_u(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_u(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_v(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_v(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_t(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_t(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_q(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_q(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_psfc(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_psfc(:,:),grid%ts_buf_size*grid%max_ts_locs)

#if (EM_CORE == 1)
   ts_buf(:,:) = grid%ts_glw(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_glw(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_gsw(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_gsw(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_hfx(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_hfx(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_lh(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_lh(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_clw(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_clw(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_rainc(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_rainc(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_rainnc(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_rainnc(:,:),grid%ts_buf_size*grid%max_ts_locs)
#endif

   ts_buf(:,:) = grid%ts_tsk(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_tsk(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_tslb(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_tslb(:,:),grid%ts_buf_size*grid%max_ts_locs)

   DEALLOCATE(ts_buf)
#endif

   IF ( wrf_dm_on_monitor() ) THEN

      iunit = get_unused_unit()
      IF ( iunit <= 0 ) THEN
         CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
      END IF

      DO i=1,grid%ntsloc_domain

         ix = grid%itsloc(i)
         iy = grid%jtsloc(i)

         OPEN(UNIT=iunit, FILE=TRIM(grid%ts_filename(i)), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

#if (EM_CORE == 1)
            WRITE(UNIT=iunit,FMT='(i2,f13.6,i5,i5,i5,1x,14(f13.5,1x))')  &
                              grid%id, grid%ts_hour(n,i),        &
                              grid%id_tsloc(i), ix, iy,          &
                              grid%ts_t(n,i),                    &
                              grid%ts_q(n,i),                    &
                              grid%ts_u(n,i),                    &
                              grid%ts_v(n,i),                    &
                              grid%ts_psfc(n,i),                 &
                              grid%ts_glw(n,i),                  &
                              grid%ts_gsw(n,i),                  &
                              grid%ts_hfx(n,i),                  &
                              grid%ts_lh(n,i),                   &
                              grid%ts_tsk(n,i),                  &
                              grid%ts_tslb(n,i),                 &
                              grid%ts_rainc(n,i),                &
                              grid%ts_rainnc(n,i),               &
                              grid%ts_clw(n,i)
#else
            WRITE(UNIT=iunit,FMT='(i2,f13.6,i5,i5,i5,1x,7(f13.5,1x))')  &
                              grid%id, grid%ts_hour(n,i),        &
                              grid%id_tsloc(i), ix, iy,          &
                              grid%ts_t(n,i),                    &
                              grid%ts_q(n,i),                    &
                              grid%ts_u(n,i),                    &
                              grid%ts_v(n,i),                    &
                              grid%ts_psfc(n,i),                 &
                              grid%ts_tsk(n,i),                  &
                              grid%ts_tslb(n,i)
#endif
         END DO

         CLOSE(UNIT=iunit)

      END DO

   END IF

   grid%next_ts_time = 1

END SUBROUTINE write_ts


#if (EM_CORE == 1)
SUBROUTINE calc_p8w(grid, ix, iy, p8w, k_start, k_end)

   USE module_domain
   USE module_model_constants

   IMPLICIT NONE

   ! Arguments
   TYPE (domain), INTENT(IN) :: grid
   INTEGER, INTENT(IN) :: ix, iy, k_start, k_end
   REAL, DIMENSION(k_start:k_end), INTENT(OUT) :: p8w

   ! Local variables
   INTEGER :: k
   REAL    :: z0, z1, z2, w1, w2 
   REAL, DIMENSION(k_start:k_end)   :: z_at_w
   REAL, DIMENSION(k_start:k_end-1) :: z


   DO k = k_start, k_end
      z_at_w(k) = (grid%phb(ix,k,iy)+grid%ph_2(ix,k,iy))/g
   END DO

   DO k = k_start, k_end-1
      z(k) = 0.5*(z_at_w(k) + z_at_w(k+1))
   END DO

   DO k = k_start+1, k_end-1
      p8w(k) = grid%fnm(k)*(grid%p(ix,k,iy)+grid%pb(ix,k,iy)) + &
               grid%fnp(k)*(grid%p(ix,k-1,iy)+grid%pb(ix,k-1,iy))
   END DO

   z0 = z_at_w(k_start)
   z1 = z(k_start)
   z2 = z(k_start+1)
   w1 = (z0 - z2)/(z1 - z2)
   w2 = 1. - w1
   p8w(k_start) = w1*(grid%p(ix,k_start,iy)+grid%pb(ix,k_start,iy)) + &
                  w2*(grid%p(ix,k_start+1,iy)+grid%pb(ix,k_start+1,iy))

   z0 = z_at_w(k_end)
   z1 = z(k_end-1)
   z2 = z(k_end-2)
   w1 = (z0 - z2)/(z1 - z2)
   w2 = 1. - w1
   p8w(k_end) = exp(w1*log(grid%p(ix,k_end-1,iy)+grid%pb(ix,k_end-1,iy)) + &
                    w2*log(grid%p(ix,k_end-2,iy)+grid%pb(ix,k_end-2,iy)))

END SUBROUTINE calc_p8w
#endif
