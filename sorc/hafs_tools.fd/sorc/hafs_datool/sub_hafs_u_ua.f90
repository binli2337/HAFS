!========================================================================================
  subroutine hafs_u_ua(u_ua, in_grid, in_file, out_file)

!-----------------------------------------------------------------------------
! HAFS DA tool - u_ua
! authors and history:
!      -- 202411, Yonghui Weng added this subroutine to update ua/va(u/v) with u/v(ua/va) values.
!
!-----------------------------------------------------------------------------
! This subroutine updates ua/va(u/v) with u/v(ua/va) values.
! u_ua:    string, "u_update_ua" or "ua_update_u"
!          u_update_ua -- update ua/va with u/v values
!          ua_update_u -- update u/v with ua/va values
! in_grid: domain grid file, grid_spec.nc or grid_mspec.nest02_2024_06_30_18.tile2.nc
! in_file: hafs restart file, usually is 20240630.180000.fv_core.res.tile1.nc or 20240630.180000.fv_core.res.nest02.tile2.nc
!-----------------------------------------------------------------------------

  use constants
  use netcdf
  use module_mpi
  use var_type

  implicit none
  character (len=*), intent(in)         :: u_ua, in_grid, in_file, out_file

  type(grid2d_info)                     :: grid
  real, allocatable, dimension(:,:)     :: cangu, sangu, cangv, sangv
  integer  :: ix, iy, iz, ndims, k, ncid, varid, i, j, rcode
  integer  :: yaxis_1, yaxis_2, xaxis_1, xaxis_2
  integer  :: yaxis_1id, yaxis_2id, xaxis_1id, xaxis_2id, uid, vid, zaxis_1id, timeid
  integer, dimension(nf90_max_var_dims) :: dims
  real, allocatable, dimension(:,:,:,:) :: dat41, dat42, ua, va, u, v
  real, allocatable, dimension(:,:)     :: dat21, dat22

!------------------------------------------------------------------------------
! 1 --- arg process and parameters

!------------------------------------------------------------------------------
! 2 --- input grid info: read from grid file grid_spec.nc
  call rd_grid_spec_data(trim(in_grid), grid)
  ix=grid%grid_xt
  iy=grid%grid_yt
  allocate( cangu(ix,iy+1),sangu(ix,iy+1),cangv(ix+1,iy),sangv(ix+1,iy) )
  call cal_uv_coeff_fv3(ix, iy, grid%grid_lat, grid%grid_lon, cangu, sangu, cangv, sangv)

!------------------------------------------------------------------------------
! 3 --- updates
  if ( trim(u_ua) == 'u_update_ua' ) then
     !---get u dimensions
     call get_var_dim(trim(in_file), 'u', ndims, dims)
     if ( ndims /= 4 ) then
        write(*,'(a)')' !!!!! warning: please check u dimensions in '//trim(in_file)
        write(*,*)ndims, dims
        stop
     endif
     iz=dims(3)
     write(*,*),'dims: ',dims(1:4)
     write(*,'(a,i,a,i,a,i,a,i,a,i)')'u dimensions: ',ix,'=',dims(1),':',iy+1,'=',dims(2),':',iz

     !---get u/v
     allocate(dat41(ix, iy+1, iz,1), dat42(ix+1, iy, iz,1))
     call get_var_data(trim(in_file), 'u', ix, iy+1, iz, 1, dat41)
     call get_var_data(trim(in_file), 'v', ix+1, iy, iz, 1, dat42)

     !---convert u,v from fv3grid to earth
     allocate(dat21(ix, iy+1), dat22(ix+1, iy))
     allocate(ua(ix,iy,iz,1), va(ix,iy,iz,1))
     do k = 1, iz
        call fv3uv2earth(ix, iy, dat41(:,:,k,1), dat42(:,:,k,1), cangu, sangu, cangv, sangv, dat21, dat22)
        !---destage: C-/D- grid to A-grid
        ua(:,:,k,1) = (dat21(:,1:iy)+dat21(:,2:iy+1))/2.0
        va(:,:,k,1) = (dat22(1:ix,:)+dat22(2:ix+1,:))/2.0
     enddo
     deallocate(dat41, dat42, dat21, dat22, cangu, sangu, cangv, sangv)

     !---update in_file
     call update_hafs_restart(trim(in_file), 'ua', ix, iy, iz, 1, ua)
     call update_hafs_restart(trim(in_file), 'va', ix, iy, iz, 1, va)
     deallocate(ua, va)

  else if ( trim(u_ua) == 'ua_update_u' ) then
     !---get ua dimensions
     call get_var_dim(trim(in_file), 'ua', ndims, dims)
     if ( ndims /= 4 ) then
        write(*,'(a)')' !!!!! warning: please check ua dimensions in '//trim(in_file)
        write(*,*)ndims, dims
        stop
     endif
     iz=dims(3)
     write(*,'(a,i,a,i,a,i,a,i,a,i)')'u dimensions: ',ix,'=',dims(1),':',iy,'=',dims(2),':',iz

     !---get ua/va
     allocate(dat41(ix, iy, iz, 1), dat42(ix, iy, iz, 1))
     call get_var_data(trim(in_file), 'ua', ix, iy, iz, 1, dat41)
     call get_var_data(trim(in_file), 'va', ix, iy, iz, 1, dat42)

     !---stage ua/va grid
     allocate(ua(ix, iy+1, iz, 1), va(ix+1, iy, iz, 1))
     ua(:,1   ,:,1)=dat41(:,1,:,1)
     ua(:,2:iy,:,1)=0.5*(dat41(:,1:iy-1,:,1)+dat41(:,2:iy,:,1))
     ua(:,iy+1,:,1)=dat41(:,iy,:,1)
     va(1   ,:,:,1)=dat42(1,:,:,1)
     va(2:ix,:,:,1)=0.5*(dat42(1:ix-1,:,:,1)+dat42(2:ix,:,:,1))
     va(ix+1,:,:,1)=dat42(ix,:,:,1)
     deallocate(dat41, dat42)

     !---convert earth wind to fv3grid wind
     allocate(u(ix, iy+1, iz,1), v(ix+1, iy, iz,1))
     do k = 1, iz
        call earthuv2fv3(ix, iy, ua(:,:,k,1), va(:,:,k,1), cangu, sangu, cangv, sangv, u(:,:,k,1), v(:,:,k,1))
     enddo
     deallocate(ua, va, cangu, sangu, cangv, sangv)

     !---check dimensions: xaxis_1=ix, xaxis_2=ix+1, yaxis_2=iy, yaxis_1=iy+1, zaxis_1=iz, zaxis_2=iz+1
     xaxis_1=-99; xaxis_2=-99; yaxis_1=-99; yaxis_2=-99
     call nccheck(nf90_open(trim(in_file), nf90_write, ncid), 'wrong in open '//trim(in_file), .true.)

     !---change dim yaxis_1 to yaxis_2
     rcode=nf90_inq_dimid(ncid, "yaxis_1", yaxis_1id)
     if ( rcode == nf90_noerr ) rcode=nf90_inquire_dimension(ncid, yaxis_1id, len=yaxis_1)
     rcode=nf90_inq_dimid(ncid, "yaxis_2", yaxis_2id)
     if ( rcode == nf90_noerr ) rcode=nf90_inquire_dimension(ncid, yaxis_2id, len=yaxis_2)
     rcode=nf90_redef(ncid)
     if ( yaxis_1 == iy .and. yaxis_2 < 1 ) then
        write(*,*)'====w1 rename dimension yaxis_1 --> yaxis_2'
        rcode=nf90_rename_dim(ncid, yaxis_1id, "yaxis_2")
        write(*,*)'====w1 rename dimension yaxis_1 --> yaxis_2 finished'
        rcode=nf90_def_dim(ncid, "yaxis_1", iy+1, yaxis_1id)
     endif

     rcode=nf90_inq_dimid(ncid, "xaxis_2", xaxis_2id)
     if ( rcode /= nf90_noerr ) rcode=nf90_def_dim(ncid, "xaxis_2", ix+1, xaxis_2id)

     !---def u/v var
     !rcode=nf90_inq_dimid(ncid, "xaxis_1", xaxis_1id)
     !rcode=nf90_inq_dimid(ncid, "yaxis_1", yaxis_1id)
     !rcode=nf90_inq_dimid(ncid, "yaxis_2", yaxis_2id)
     !rcode=nf90_inq_dimid(ncid, "zaxis_1", zaxis_1id)
     !rcode=nf90_inq_dimid(ncid, "Time", timeid)
     !rcode=nf90_def_var(ncid, 'u', nf90_double, (/xaxis_1id,yaxis_1id,zaxis_1id,timeid/), uid)
     !rcode=nf90_def_var(ncid, 'v', nf90_double, (/xaxis_2id,yaxis_2id,zaxis_1id,timeid/), vid)
     rcode=nf90_enddef(ncid)

     !---output u/v
     !rcode=nf90_put_var(ncid, uid, dble(u))
     !rcode=nf90_put_var(ncid, vid, dble(v))
     rcode=nf90_close(ncid)

     !      rcode=nf90_inq_dimid(ncid, "yaxis_1", yaxis_1id)
     !      rcode=nf90_rename_dim(ncid, yaxis_1id, "yaxis_2")
           !rcode=nf90_inq_varid(ncid, 'yaxis_1', varid)
           !rcode=nf90_rename_var(ncid, varid, "yaxis_2")
           !rcode=nf90_rename_att(ncid, varid, "long_name", "yaxis_2")
     !   endif
     !   call nccheck(nf90_enddef(ncid), 'wrong in nf90_enddef '//trim(in_file), .true.)
     !endif
     !call nccheck(nf90_close(ncid), 'wrong in close '//trim(in_file), .true.)

     !call nccheck(nf90_open(trim(in_file), nf90_write, ncid), 'wrong in open '//trim(in_file), .true.)
     !rcode=nf90_inq_varid(ncid, 'yaxis_1', varid)
     !rcode=nf90_rename_var(ncid, varid, "yaxis_2")
     !rcode=nf90_rename_att(ncid, varid, "long_name", "yaxis_2")
     !call nccheck(nf90_close(ncid), 'wrong in close '//trim(in_file), .true.)

     !call nccheck(nf90_inq_dimid(ncid, "xaxis_1", xaxis_1id), 'wrong in nf90_inq_dimid xaxis_1', .false.)
     !if ( xaxis_1id > 0 ) call nccheck(nf90_inquire_dimension(ncid, xaxis_1id, len=xaxis_1), 'wrong in nf90_inquire_dimension xaxis_1', .false.)
     !call nccheck(nf90_inq_dimid(ncid, "xaxis_2", xaxis_2id), 'wrong in nf90_inq_dimid xaxis_2', .false.)
     !if ( xaxis_2id > 0 ) call nccheck(nf90_inquire_dimension(ncid, xaxis_2id, len=xaxis_2), 'wrong in nf90_inquire_dimension xaxis_2', .false.)
     !if ( xaxis_1/=ix .or. xaxis_2/=ix+1 .or. yaxis_2/=iy .or. yaxis_1/=iy+1 ) then
     !      write(*,*)'====w5', iy, yaxis_1, yaxis_2
     !      call nccheck(nf90_def_dim(ncid, "yaxis_1", iy+1, i), 'wrong in nf90_def_dim yaxis_1', .true.)
     !      write(*,*)'====w3 redef dim yaxis_1'
     !   endif
     !   !write(*,*)'====w6', ix, xaxis_1, xaxis_2
     !   !if ( xaxis_2 < 1 ) then
     !   !   call nccheck(nf90_def_dim(ncid, "xaxis_2", ix+1, xaxis_2id), 'wrong in nf90_inq_dimid xaxis_2', .true.)
     !   !   write(*,*)'====w4 def dim xaxis_2'
     !   !endif
     !   call nccheck(nf90_enddef(ncid), 'wrong in nf90_enddef '//trim(in_file), .true.)
     !endif
     !call nccheck(nf90_close(ncid), 'wrong in close '//trim(in_file), .true.)
     !write(*,*)'====w1 finished redef'

     !---update in_file
     !call update_hafs_restart(trim(in_file), 'u', ix, iy+1, iz, 1, u)
     !call update_hafs_restart(trim(in_file), 'v', ix+1, iy, iz, 1, v)
     !call write_nc_real(trim(in_file),'u', ix, iy+1, iz, 1, 'xaxis_1', 'yaxis_1', 'zaxis_1', 'Time', u, 'm-1s', 'u')
     !call write_nc_real(trim(in_file),'v', ix+1, iy, iz, 1, 'xaxis_2', 'yaxis_2', 'zaxis_1', 'Time', v, 'm-1s', 'v')
     deallocate(u, v)

  endif

  return
  end subroutine hafs_u_ua

!========================================================================================
