MODULE mo_lagrangian_tracker
  USE grid, ONLY : deltax, deltay, deltaz, xt,yt,zt
  USE util, ONLY : arr_resize, closest
  IMPLICIT NONE

  CONTAINS
  
  SUBROUTINE lagrangian_tracker(ix,iy,iz,t,np,t_in,t_out,t_trac,scS,start_time,end_time,x,y,z)
    
    CHARACTER(len=50), PARAMETER :: name = "lagrangian_tracker"
    
    REAl, INTENT(in) :: scS
    REAl, INTENT(in) :: start_time  
    REAl, INTENT(inout) :: end_time
    REAl, INTENT(inout) :: t_trac
    REAl, INTENT(in) :: x(:), y(:), z(:)
    INTEGER, INTENT(inout) :: np
    
    REAl, ALLOCATABLE, INTENT(out) ::  t(:),t_in(:),t_out(:)
    INTEGER, ALLOCATABLE, INTENT(out) ::  ix(:), iy(:), iz(:)
    
    INTEGER :: i, j, k, io, ix1, iy1, iz1
    INTEGER :: np_trac   
    REAl    :: deltax2, deltay2, deltaz2, d, d_x, d_y, d_z, dx, dy, dz, dx2x1, dy2y1, dz2z1, &
         x_trac, y_trac, z_trac, vx, vy, vz, xlim1, xlim2, ylim1, ylim2, x1, x2, y1, y2, z1, z2  
    REAl    :: t_trac_old, sub_dt, dt, dtx, dty, dtz, min_t, dt2(2) = -999., dt4(4) = -999.
    REAl    :: eps = 1e-6
    LOGICAL :: nexpout = .FALSE.
    
    deltax2 = deltax/2
    deltay2 = deltay/2
    deltaz2 = deltaz/2
    np_trac = 2
    ALLOCATE(ix(np_trac),iy(np_trac),iz(np_trac),t(np_trac))
    ALLOCATE(t_in(1),t_out(1))
    ix = -999
    iy = -999
    iz = -999
    t  = -999.
    t_in = -999.
    t_out = t_in - 999.
    
    xlim1 = MINVAL(xt) + 3*deltax2
    xlim2 = MAXVAL(xt) - 3*deltax2
    ylim1 = MINVAL(yt) + 3*deltay2
    ylim2 = MAXVAL(yt) - 3*deltay2
        
    x_trac = xlim1 - 999.
    y_trac = ylim1 - 999.
    z_trac = -999.
    
    j  = 1
    io = 0
    DO i = 1,np-1
       t_trac = start_time
       d_x = x(i+1) - x(i)
       d_y = y(i+1) - y(i)
       d_z = z(i+1) - z(i)
       
       IF ( (d_x > eps) .AND. (ABS(d_y) > eps) ) THEN
          dz2z1 = 0.
          IF (d_y >  eps) THEN
             x1 = (ylim1-y(i))*(d_x/d_y) + x(i)
             x1 = MAX(x1,MAX(x(i),xlim1))
             y1 = (x1-x(i))*(d_y/d_x) + y(i)
             x2 = (ylim2-y(i))*(d_x/d_y) + x(i)
          END IF
          IF (d_y < -eps) THEN
             x1 = (ylim2-y(i))*(d_x/d_y) + x(i)
             x1 = MAX(x1,MAX(x(i),xlim1))
             y1 = (x1-x(i))*(d_y/d_x) + y(i)
             x2 = (ylim1-y(i))*(d_x/d_y) + x(i)
          END IF
          x2 = MIN(x2,MIN(x(i+1),xlim2))
          y2 = (x2-x(i))*(d_y/d_x) + y(i)
          dx2x1 = x2 - x1
          dy2y1 = y2 - y1
          
       ELSE IF ( (d_x < -eps) .AND. (ABS(d_y) > eps) )  THEN
          dz2z1 = 0.
          IF (d_y >  eps) THEN
             x1 = (ylim1-y(i))*(d_x/d_y) + x(i)
             x1 = MIN(x1,MIN(x(i),xlim2))
             y1 = (x1-x(i))*(d_y/d_x) + y(i)
             y1 = (x1-x(i))*(d_y/d_x) + y(i)
             x2 = (ylim2-y(i))*(d_x/d_y) + x(i)
          END IF
          IF (d_y < -eps) THEN
             x1 = (ylim2-y(i))*(d_x/d_y) + x(i)
             x1 = MIN(x1,MIN(x(i),xlim2))
             y1 = (x1-x(i))*(d_y/d_x) + y(i)
             x2 = (ylim1-y(i))*(d_x/d_y) + x(i)
          END IF
          x2 = MAX(x2,MAX(xlim1,x(i+1)))
          y2 = (x2-x(i))*(d_y/d_x) + y(i)
          dx2x1 = x2 - x1
          dy2y1 = y2 - y1
          
       ELSE IF ( (ABS(d_x) > eps) .AND. (ABS(d_y) <= eps) )  THEN
          dz2z1 = 0.
          IF (d_x >  eps) THEN
             x1 = MAX(x(i),xlim1)
             x2 = MIN(x(i+1),xlim2)
             dx2x1 = x2 - x1
          END IF
          IF (d_x < -eps) THEN
             x1 = MIN(x(i),xlim2)
             x2 = MAX(x(i+1),xlim1)
             dx2x1 = x2 - x1
          END IF
          IF ((y(i) > ylim1).AND.(y(i) < ylim2)) THEN
             dy2y1 = 0.
             d_y = 0.
             y1 = y(i)
             y2 = y(i)
          ELSE
             dy2y1 = -1.
             d_y = 1.
          END IF
          
       ELSE IF ( (ABS(d_x) <= eps) .AND. (ABS(d_y) >  eps) )  THEN
          dz2z1 = 0.
          IF (d_y >  eps) THEN
             y1 = MAX(y(i),ylim1)
             y2 = MIN(y(i+1),ylim2)
             dy2y1 = y2 - y1 
          END IF
          IF (d_y < -eps) THEN
             y1 = MIN(y(i),ylim2)
             y2 = MAX(y(i+1),ylim1)
             dy2y1 = y2 - y1 
          END IF
          IF ((x(i) > xlim1).AND.(x(i) < xlim2)) THEN
             dx2x1 = 0.
             d_x = 0.
             x1 = x(i)
             x2 = x(i)
          ELSE
             dx2x1 = -1.
             d_x = 1.
          END IF
          
       ELSE IF ( (ABS(d_x) <= eps) .AND. (ABS(d_y) <= eps) )  THEN
          z1 = z(i)
          z2 = z(i+1)
          dz2z1 = z2 - z1
          IF ( ((x(i) > xlim1).AND.(x(i) < xlim2)).AND.((y(i) > ylim1).AND.(y(i) < ylim2)) ) THEN
             dx2x1 = 0.
             d_x = 0.
             x1 = x(i)
             x2 = x(i)
             dy2y1 = 0.
             d_y = 0.
             y1 = y(i)
             y2 = y(i)
          ELSE
             dx2x1 = -1.
             d_x = 1.
             dy2y1 = -1.
             d_y = 1.
          END IF
          
       END IF
       
       IF ((d_x*dx2x1 >= 0.).AND.(d_y*dy2y1 >= 0.).AND.(ANY([abs(dx2x1),abs(dy2y1),abs(dz2z1)] > eps))) THEN
          
          IF ((ABS(d_x) > eps) .AND. (ABS(d_y) < eps)) THEN
             z1 = (d_z/d_x) * (x1-x(i)) + z(i)
             z2 = (d_z/d_x) * (x2-x(i)) + z(i)
             dz2z1 = z2 - z1 
          ELSE IF ((ABS(d_x) < eps) .AND. (ABS(d_y) > eps)) THEN
             z1 = (d_z/d_y) * (y1-y(i)) + z(i)
             z2 = (d_z/d_y) * (y2-y(i)) + z(i)
             dz2z1 = z2 - z1
          ELSE IF ((ABS(d_x) > eps) .AND. (ABS(d_y) > eps)) THEN
             z1 = (d_z/d_y) * (y1-y(i)) + z(i)
             z2 = (d_z/d_y) * (y2-y(i)) + z(i)
             dz2z1 = z2 - z1
          END IF
          
          IF ((j == 1) .OR. ( ( x1 < (x_trac - eps) ).OR.( x1 > (x_trac + eps) ) )) THEN  
             ix1 = closest(xt,x1)
             ix(j) = ix1
          END IF
          
          IF ((j == 1) .OR. ( (y1 < (y_trac - eps) ).OR.( y1 > (y_trac + eps) ) )) THEN    
             iy1 = closest(yt,y1)
             iy(j) = iy1
          END IF
          
          IF ((j == 1) .OR. ( (z1 < (z_trac - eps) ).OR.( z1 > (z_trac + eps) ) )) THEN  
             iz1 = closest(zt,z1)
             iz(j) = iz1
          END IF
          
          nexpout = ((x(i+1) < xlim1).OR.(x(i+1) > xlim2).OR.(y(i+1) < ylim1).OR.(y(i+1) > ylim2))
          
          IF (i == 1) THEN
             t_trac = t_trac + sqrt((x1-x(i))**2 + (y1-y(i))**2 + (z1-z(i))**2)/scS
             t_trac = t_trac 
             t(j) = t_trac
             io = io + 1
             t_in(io)  = t_trac
          END IF
          
          IF (i > 1) THEN
             IF ((x_trac >= xlim1 + eps ).AND.(x_trac <= xlim2 - eps ) &
                  .AND.(y_trac >= ylim1 + eps ).AND.(y_trac<= ylim2 - eps )) THEN
                
                t_trac = t_trac_old
             ELSE
                t_trac = t_trac + sqrt((x1-x(i))**2 + (y1-y(i))**2 + (z1-z(i))**2)/scS
                IF (i > 2) THEN
                   DO k = 1,i-2
                      d_x = x(k+1) - x(k)
                      d_y = y(k+1) - y(k)
                      d_z = z(k+1) - z(k)
                      d  = sqrt(d_x**2 + d_y**2 + d_z**2)
                      dt = d/scS
                      vx = d_x/dt
                      vy = d_y/dt
                      vz = d_z/dt
                      t_trac = t_trac + dt
                      t(j)   = t_trac
                   END DO
                END IF
                d_x = x(i) - x(i-1)
                d_y = y(i) - y(i-1)
                d_z = z(i) - z(i-1)
                d  = sqrt(d_x**2 + d_y**2 + d_z**2)
                dt = d/scS
                vx = d_x/dt
                vy = d_y/dt
                
                IF ((x(i) >= xlim1 + eps ).AND.(x(i) <= xlim2 - eps ) &
                     .AND.(y(i) >= ylim1 + eps).AND.(y(i) <= ylim2 - eps)) THEN
                   
                   IF  (d_x > eps) THEN
                      dx  = xt(ix(j)) - deltax2 - x(i-1) 
                      dtx = dx/vx
                   ELSE IF (d_x < -eps) THEN
                      dx  = xt(ix(j)) + deltax2 - x(i-1) 
                      dtx = dx/vx
                   ELSE
                      dtx = -999.
                   END IF
                   
                   IF ( (d_y > eps) .AND. (y(i) >= y1) ) THEN
                      dy  = yt(iy(j)) - deltay2 - y(i-1) 
                      dty = dy/vy
                   ELSE IF ( (d_y < -eps) .AND. (y(i) <= y1) ) THEN
                      dy  = yt(iy(j)) + deltay2 - y(i-1) 
                      dty = dy/vy
                   ELSE 
                      dty = -999.
                   END IF
                   
                   IF ((ABS(dx) > eps) .OR. (ABS(dy) > eps)) THEN
                      dt2 = [dtx,dty]
                      min_t = MAXVAL(dt2, MASK = dt2 > 0.) 
                      t_trac = t_trac + min_t
                      t(j)   = t_trac
                   ELSE
                      t(j)   = t_trac
                   END IF
                   
                ELSE
                   t_trac = t_trac + dt
                   t(j)   = t_trac
                END IF
                
                IF (nexpout.OR.(x(i) < xlim1).OR.(x(i) > xlim2).OR.(y(i) < ylim1).OR.(y(i) > ylim2)) THEN
                   io = io + 1
                   CALL arr_resize(t_in,io)
                   CALL arr_resize(t_out,io)
                   t_in(io)  = t_trac
                END IF
                
             END IF
             
          END IF
          
          d  = sqrt(dx2x1**2 + dy2y1**2 + dz2z1**2)
          dt = d/scS
          vx = dx2x1/dt
          vy = dy2y1/dt
          vz = dz2z1/dt
          sub_dt = 0.
          x_trac = x1
          y_trac = y1
          z_trac = z1
          end_time = t_trac + dt
          
          DO WHILE (sub_dt < dt)
             
             IF (dx2x1 > eps)  THEN
                dx  = xt(ix(j)) + deltax2 - x_trac
                dtx = dx/vx
             ELSE IF (dx2x1 < -eps) THEN
                dx  = xt(ix(j)) - deltax2 - x_trac
                dtx = dx/vx
             ELSE
                dtx = -999.
                dx = 0.
                ix(j+1) = ix(j)
             END IF
             
             IF (dy2y1 > eps)  THEN
                dy  = yt(iy(j)) + deltay2 - y_trac
                dty = dy/vy
             ELSE IF (dy2y1 < -eps) THEN
                dy  = yt(iy(j)) - deltay2 - y_trac
                dty = dy/vy
             ELSE
                dty = -999.
                dy = 0. 
                iy(j+1) = iy(j)
             END IF
             
             IF (dz2z1 > eps)  THEN
                dz  = zt(iz(j)) + deltaz2 - z_trac
                dtz = dz/vz
             ELSE IF (dz2z1 < -eps) THEN
                dz  = zt(iz(j)) - deltaz2 - z_trac
                dtz = dz/vz
             ELSE
                dtz = -999.
                dz = 0. 
                iz(j+1) = iz(j) 
             END IF
             
             IF (d > eps)  THEN
                dt4 = [dtx,dty,dtz,(dt-sub_dt)]
                min_t = MINVAL(dt4, MASK = dt4 >= 0.)
                
                IF (min_t == dtx) THEN
                   IF (vx > 0)  ix(j+1) = ix(j) + 1
                   IF (vx < 0)  ix(j+1) = ix(j) - 1
                   IF (ABS(dtx-dty) < eps) THEN
                      IF(vy > 0) iy(j+1) = iy(j) + 1
                      IF(vy < 0) iy(j+1) = iy(j) - 1
                   ELSE
                      iy(j+1) = iy(j)
                   END IF
                   
                   iz(j+1) = iz(j)
                   t_trac = t_trac + min_t
                   t(j+1) = t_trac
                   t_trac_old = t_trac
                   
                ELSE IF (min_t == dty) THEN
                   IF (ABS(dtx-dty) < eps) THEN
                      IF(vx > 0) ix(j+1) = ix(j) + 1
                      IF(vx < 0) ix(j+1) = ix(j) - 1
                   ELSE
                      ix(j+1) = ix(j)
                   END IF
                   IF (vy > 0)  iy(j+1) = iy(j) + 1
                   IF (vy < 0)  iy(j+1) = iy(j) - 1
                   iz(j+1) = iz(j)
                   t_trac = t_trac + min_t
                   t(j+1) = t_trac
                   t_trac_old = t_trac
                   
                ELSE IF (min_t == dtz) THEN
                   ix(j+1) = ix(j) 
                   iy(j+1) = iy(j)
                   IF (vz > 0)  iz(j+1) = iz(j) + 1
                   IF (vz < 0)  iz(j+1) = iz(j) - 1
                   t_trac = t_trac + min_t
                   t(j+1) = t_trac
                   t_trac_old = t_trac           
                   
                ELSE IF (min_t == (dt-sub_dt)) THEN
                   t_trac = t_trac + min_t
                   t(j+1) = t_trac
                   t_trac_old = t_trac
                   IF ((ABS(dtx-min_t) < eps).AND.(ABS(dty-min_t) < eps)) THEN
                      IF(vx > 0) ix(j+1) = ix(j) + 1
                      IF(vx < 0) ix(j+1) = ix(j) - 1
                      IF(vy > 0) iy(j+1) = iy(j) + 1
                      IF(vy < 0) iy(j+1) = iy(j) - 1
                   ELSE IF (ABS(dtx-min_t) < eps) THEN
                      IF(vx > 0) ix(j+1) = ix(j) + 1
                      IF(vx < 0) ix(j+1) = ix(j) - 1
                   ELSE IF (ABS(dty-min_t) < eps) THEN
                      IF(vy > 0) iy(j+1) = iy(j) + 1
                      IF(vy < 0) iy(j+1) = iy(j) - 1
                   ELSE
                      j = j - 1
                   END IF
                   
                END IF
                x_trac = x_trac + vx * min_t
                y_trac = y_trac + vy * min_t
                z_trac = z_trac + vz * min_t
                sub_dt = sub_dt + min_t    
             END IF
             
             IF ( j+1 == np_trac) THEN
                np_trac = (j+1) * 2
                CALL arr_resize(ix,np_trac)
                CALL arr_resize(iy,np_trac)
                CALL arr_resize(iz,np_trac)
                CALL arr_resize(t,np_trac)
             END IF
             
             IF (t(j+1) == t(j)) THEN  
                ix(j) = ix(j+1)
                iy(j) = iy(j+1)
                iz(j) = iz(j+1)
                j = j - 1
             END IF
             
             j = j + 1
             
          END DO
          
          IF ( ((io >= 1).AND.nexpout).OR.(i+1 == np)) THEN
             t_out(io) = end_time
          END IF
          
       END IF
       
    END DO
    
    t_trac = t_in(1)
    np = j
    IF (j > 1) THEN
       CALL arr_resize(ix,np)
       CALL arr_resize(iy,np)
       CALL arr_resize(iz,np)
       CALL arr_resize(t,np+1)
       t(np+1) = t_out(io)
    END IF
    
  END SUBROUTINE lagrangian_tracker
  
  
END MODULE mo_lagrangian_tracker
