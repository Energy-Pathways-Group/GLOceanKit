subroutine write_results
 use io_params 
 use etc,           only: istep,istart 
 use particles,     only: nparticles,pinc=>particle_write_inc
 use mpi_params,    only: myid,comm,ierr
  
 implicit none
 integer               :: fid
 integer,save          :: jstart(maxsets)
 logical,save          :: first_entry=.TRUE.
 
 if( first_entry ) then
  jstart(:)=istart      ! generally, all output sets start up at istep=istart
  first_entry=.false.
 endif
  
 if( istep >= 0 .and. nparticles > 0 .and. mod(istep,pinc)==0 ) then 
  call write_particle_data
 endif
 
 !-------------------------------------------------------
 ! loop through all the user-specified file sets...
 !-------------------------------------------------------
 do fid=1,num_file_sets 
 
  if( istep >= jstart(fid) ) then
  
   if( mod( (istep-jstart(fid)),nsteps(fid) ) == 0    &
       .or. istep == jstart(fid) ) then
    
    !-------------------------------------------------------
    ! initialize netcdf file if necessary
    !-------------------------------------------------------
    if( trim(mode(fid)) == 'new' .or. istep == jstart(fid)) then      
     call init_netcdf(fid) 
    endif 
   
    !-------------------------------------------------------
    ! file already initialized, just write current data
    !-------------------------------------------------------
    call write_netcdf(fid)
   
   endif 
  
  endif
  
 enddo   ! end loop through file sets

 call mpi_barrier(comm,ierr)
 return
end subroutine write_results



  
  
subroutine init_netcdf(fid)
 use io_params
 use decomposition_params
 use dimensional_scales,     only: length_scale,scalar_scale
 use methods_params,         only: do_second_scalar
 use dependent_variables
 use independent_variables
 use etc,                    only: istep
 use mpi_params
      
 implicit none
 include 'netcdf.inc'

 character(len=80)              :: filename
 character(len=3)               :: dimstr
 character(len=80),save         :: topdir='output/'
 character(len=80)              :: s1_name,s2_name
 character(len=80)              :: s1_units,s2_units
 integer                        :: i,j,k
 integer                        :: fid,ncid,rcode
 integer                        :: xVarID,iid
 integer                        :: yVarID,jid
 integer                        :: zVarID,kid
 integer                        :: tVarID,timeid
 integer                        :: All_nD_VarSIZE(4),ndims
 integer                        :: uVarID,vVarID,wVarID 
 integer                        :: s1VarID,s2VarID
 integer                        :: s1_barVarID,s2_barVarID
 integer                        :: divustarVarID,phiVarID,pdVarID
 integer                        :: s1_name_len,s2_name_len
 integer                        :: s1_units_len,s2_units_len
 integer                        :: npts,start1D(1),count1D(1)
 integer                        :: counter,offset
 real(kind=4),allocatable,save  :: scratch(:)
 logical,save                   :: first_entry=.TRUE.

 if( first_entry ) then
  npts = maxval((/nx,ny,nz/))
  allocate( scratch(npts)  )
  first_entry=.FALSE.
 endif
 
 call process_dimensions(fid,dimstr)  !! figure out if file should be 1D,2D,3D etc
 
 call construct_filename(fid,dimstr,topdir,istep,fullname(fid))
 
 !---------------------------------------------------------------
 ! specify the start and count arrays used in netcdf write calls
 ! set the logical variable do_write(fid), this routine will assume
 ! that there is data to write for myid and detect if this is
 ! really the case. If not, it will set do_write(fid) to .FALSE.
 !---------------------------------------------------------------
 call set_local_indices_count(fid)
 time_counter(fid)=1    ! always 1 for initialization
 if( .NOT. do_write(fid) ) goto 999
 
 
 !--------------------------------------------------------------- 
 !  Open a netcdf file
 !---------------------------------------------------------------
 rcode=NF_CREATE(trim(fullname(fid)),NF_NOCLOBBER,ncid)
 if(rcode.ne.NF_NOERR) then
  write(0,*) '... ERROR OPENING NETCDF FILE: init_netcdf ',trim(fullname(fid))
  write(0,*) '... myid, rcode ',myid, rcode
  stop
 endif
  

!--------------------------------------------------------------- 
! Define (and order) the spatial dimensions
!--------------------------------------------------------------- 
 rcode=NF_DEF_DIM(ncid,'idimension',count(dimid_x(fid),fid),iid)
 if (rcode.ne.NF_NOERR) write(0,*) myid,  &
         ': NetCDF Error: NF_DEF_DIM: idimension', rcode
  All_nD_VarSIZE(dimid_x(fid))=iid

  rcode=NF_DEF_DIM(ncid,'jdimension',count(dimid_y(fid),fid),jid)
  if (rcode.ne.NF_NOERR) write(0,*) myid,  &
          ': NetCDF Error: NF_DEF_DIM: jdimension', rcode
  All_nD_VarSIZE(dimid_y(fid))=jid

  rcode=NF_DEF_DIM(ncid,'kdimension',count(dimid_z(fid),fid),kid)
  if (rcode.ne.NF_NOERR) write(0,*) myid,   &
          ': NetCDF Error: NF_DEF_DIM: kdimension', rcode
  All_nD_VarSIZE(dimid_z(fid))=kid
   
  rcode=NF_DEF_DIM(ncid,'timedimension',NF_UNLIMITED,timeid)
  if (rcode.ne.NF_NOERR) write(0,*) myid,   &
          ': NetCDF Error: NF_DEF_DIM: timedimension', rcode
  All_nD_VarSIZE(dimid_t(fid))=timeid
  ndims = nspace(fid) + 1
  

!---------------------------------------------------------------
!  Define x,y,z grid position  and time variables
!---------------------------------------------------------------

!!**X*****
  rcode=NF_DEF_VAR(ncid,'x',NF_FLOAT,1,iid,xVarID)
  if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_DEF_VAR: x', rcode
  rcode=NF_PUT_ATT_TEXT(ncid,xVarID,'long_name',1,'x')
  if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_PUT_ATT_TEXT: x', rcode
  rcode=NF_PUT_ATT_TEXT(ncid,xVarID,'units',1,'m')
  if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_PUT_ATT_TEXT: x', rcode

!!**Y*****
  rcode=NF_DEF_VAR(ncid,'y',NF_FLOAT,1,jid,yVarID)
  if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_DEF_VAR: y', rcode
  rcode=NF_PUT_ATT_TEXT(ncid,yVarID,'long_name',1,'y')
  if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_PUT_ATT_TEXT: y', rcode
  rcode=NF_PUT_ATT_TEXT(ncid,yVarID,'units',1,'m')
  if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_PUT_ATT_TEXT: y', rcode
           
!!**Z*****
  rcode=NF_DEF_VAR(ncid,'z',NF_FLOAT,1,kid,zVarID)
  if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_DEF_VAR: z', rcode
  rcode=NF_PUT_ATT_TEXT(ncid,zVarID,'long_name',1,'z')
  if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_PUT_ATT_TEXT: z', rcode
  rcode=NF_PUT_ATT_TEXT(ncid,zVarID,'units',1,'m')
  if (rcode.ne.NF_NOERR) print *,myid,   &
       ': NetCDF Error: NF_PUT_ATT_TEXT: z', rcode

!!**time*****
  rcode=NF_DEF_VAR(ncid,'time',NF_FLOAT,1,timeid,tVarID)
  if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_DEF_VAR: time', rcode
  rcode=NF_PUT_ATT_TEXT(ncid,tVarID,'long_name',4,'time')
  if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_PUT_ATT_TEXT: time', rcode
  rcode=NF_PUT_ATT_TEXT(ncid,tVarID,'units',7,'seconds')
  if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_PUT_ATT_TEXT: time', rcode     

  
!-----------------------------------------------
! Set the appropriate labels for s1 and s2.
!-----------------------------------------------
  if(trim(scalar_kind(1)) == 't') then
    s1_name = 'Temperature'
    s1_name_len = 11
    s1_units = 'deg C'
    s1_units_len = 5
  endif
  if(trim(scalar_kind(2)) == 't') then
    s2_name = 'Temperature'
    s2_name_len = 11
    s2_units = 'deg C'
    s2_units_len = 5
  endif

  if(trim(scalar_kind(1)) == 's') then
    s1_name = 'Salinity'
    s1_name_len = 8
    s1_units = 'psu'
    s1_units_len = 3
  endif
  if(trim(scalar_kind(2)) == 's') then
    s2_name = 'Salinity'
    s2_name_len = 8
    s2_units = 'psu'
    s2_units_len = 3
  endif
       
  if(trim(scalar_kind(1)) == 'p') then
    s1_name = 'Passive Tracer'
    s1_name_len = 14
    s1_units = 'Concentration'
    s1_units_len = 13
  endif
  if(trim(scalar_kind(2)) == 'p') then
    s2_name = 'Passive Tracer'
    s2_name_len = 14
    s2_units = 'Concentration'
    s2_units_len = 13
  endif
             
  if(trim(scalar_kind(1)) == 'r') then
    s1_name = 'Density'
    s1_name_len = 7
    s1_units = 'kg/m3'
    s1_units_len = 5
  endif
  if(trim(scalar_kind(2)) == 'r') then
    s1_name = 'Density'
    s1_name_len = 7
    s1_units = 'kg/m3'
    s1_units_len = 5
  endif
  
    

!!
! scalar (1)
!!   
    if( variable_key(4,fid) /= 0 ) then    
      rcode=NF_DEF_VAR(ncid,'s1_bar',NF_FLOAT,1,kid,s1_barVarID)
      if (rcode.ne.NF_NOERR) print *,myid,   &
          ': NetCDF Error: NF_DEF_VAR: s1_bar', rcode       
      rcode=NF_DEF_VAR(ncid,'s1',NF_FLOAT,ndims,All_nD_VarSIZE(1:ndims),s1VarID)
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: s1', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,s1VarID,'long_name',s1_name_len,s1_name)
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: s1', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,s1VarID,'units',s1_units_len,s1_units)
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: s1', rcode
    endif

!!
!scalar (2)
!!
    if( do_second_scalar .and. variable_key(5,fid) /= 0 ) then
      rcode=NF_DEF_VAR(ncid,'s2_bar',NF_FLOAT,1,kid,s2_barVarID)
      if (rcode.ne.NF_NOERR) print *,myid,   &
          ': NetCDF Error: NF_DEF_VAR: s2_bar', rcode
      rcode=NF_DEF_VAR(ncid,'s2',NF_FLOAT,ndims,All_nD_VarSIZE,s2VarID)
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: s2', rcode           
      rcode=NF_PUT_ATT_TEXT(ncid,s2VarID,'long_name',s2_name_len,s2_name)                           
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: s2', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,s2VarID,'units',s2_units_len,s2_units)
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: s2', rcode
     endif

!!
!x velocity = u
!!
    if( variable_key(1,fid) == 1 ) then
      rcode=NF_DEF_VAR(ncid,'u',NF_FLOAT,ndims,All_nD_VarSIZE,uVarID)
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: u', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,uVarID,'long_name',1,'u')
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: u', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,uVarID,'units',3,'m/s')
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: u', rcode
    endif
!!
!y velocity = v
!! 
    if( variable_key(2,fid) == 1 ) then
      rcode=NF_DEF_VAR(ncid,'v',NF_FLOAT,ndims,All_nD_VarSIZE,vVarID)
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: v', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,vVarID,'long_name',1,'v')
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: v', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,vVarID,'units',3,'m/s')
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: v', rcode
    endif
!!
!z velocity = w
!!
    if( variable_key(3,fid) == 1 ) then
      rcode=NF_DEF_VAR(ncid,'w',NF_FLOAT,ndims,All_nD_VarSIZE,wVarID)
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: w', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,wVarID,'long_name',1,'w')
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: w', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,wVarID,'units',3,'m/s')
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: w', rcode
    endif
        
!!
!divustar = divustar
!! 
    if( variable_key(6,fid) == 1 ) then
      rcode=NF_DEF_VAR(ncid,'divustar',NF_FLOAT,ndims,All_nD_VarSIZE,divustarVarID)
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: divustar', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,divustarVarID,'long_name',9,'div ustar')                            
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: divustar', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,divustarVarID,'units',4,'1/s')
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: divustar', rcode
    endif
    
!!
!pressure variable = phi
!! 
    if( variable_key(7,fid) == 1 ) then
      rcode=NF_DEF_VAR(ncid,'phi',NF_FLOAT,ndims,All_nD_VarSIZE,phiVarID)
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: phi', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,phiVarID,'long_name',8,'pressure')                            
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: phi', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,phiVarID,'units',6,'kg/ms2')
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: phi', rcode
    endif
    
!!
!perterbation density = pd
!! 
    if( variable_key(8,fid) == 1 ) then
      rcode=NF_DEF_VAR(ncid,'pdVar',NF_FLOAT,ndims,All_nD_VarSIZE,pdVarID)
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: pd', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,pdVarID,'long_name',20,'perturbation density')                            
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: pd', rcode
      rcode=NF_PUT_ATT_TEXT(ncid,pdVarID,'units',5,'kg/m3')
       if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: pd', rcode
    endif
    
!------------------------------------------------------------------------
!End define mode
!------------------------------------------------------------------------
    rcode=NF_ENDDEF(ncid)
    if(rcode.ne.NF_NOERR) print *,myid,'ERROR  LEAVING DEFINE MODE'
!------------------------------------------------------------------------
!End define mode
!------------------------------------------------------------------------

!    store x values in temp array
      counter = 1
      offset = global_x_indices(START,YBLOCK,myid) - 1
      do i = my_x0(fid)+offset,my_x1(fid)+offset,my_xinc(fid)
       scratch(counter) = x(i)*length_scale ! [m]
       counter = counter + 1
      enddo

!    Put dimension data in netcdf file
      start1D(1) = 1
      count1D(1) = my_nx(fid)
      rcode=NF_PUT_VARA_REAL(ncid,xVarID,start1D,count1D,scratch)
      if(rcode.ne.NF_NOERR) print *,myid,'NetCDF ERROR xVarID: ',rcode


!    store y values in temp array
      counter = 1
      offset = global_y_indices(START,YBLOCK,myid) - 1
      do j = my_y0(fid)+offset,my_y1(fid)+offset,my_yinc(fid)
       scratch(counter) = y(j)*length_scale ! [m]
       counter = counter + 1
      enddo

!    Put dimension data in netcdf file
      start1D(1) = 1
      count1D(1) = my_ny(fid)
      rcode=NF_PUT_VARA_REAL(ncid,yVarID,start1D,count1D,scratch)
      if(rcode.ne.NF_NOERR) print *,myid,'NetCDF ERROR yVarID: ',rcode,filename,count1D(1)

    

!    store z values in temp array
      counter = 1
      offset = global_z_indices(START,YBLOCK,myid) - 1
      do k = my_z0(fid)+offset,my_z1(fid)+offset,my_zinc(fid)
       scratch(counter) = z(k)*length_scale ! [m]
       counter = counter + 1
      enddo

!    Put dimension data in netcdf file
      start1D(1) = 1
      count1D(1) = my_nz(fid)
      rcode=NF_PUT_VARA_REAL(ncid,zVarID,start1D,count1D,scratch)
      if(rcode.ne.NF_NOERR) print *,myid,'NetCDF ERROR zVarID: ',rcode
      

!     Write s1_bar if s1 itself is to be written     
      if( variable_key(4,fid) /= 0 ) then
       counter = 1
       offset = global_z_indices(START,YBLOCK,myid) - 1
       do k = my_z0(fid)+offset,my_z1(fid)+offset,my_zinc(fid)
        scratch(counter) = s1_bar(k,1)*scalar_scale(1)
        counter = counter + 1
       enddo
       rcode=NF_PUT_VARA_REAL(ncid,s1_barVarID,start1D,count1D,scratch)
       if(rcode.ne.NF_NOERR) print *,myid,'NetCDF ERROR s1_barVarID: ',rcode
      endif
      
!     Write s2_bar if s2 itself is to be written     
      if( do_second_scalar .and. variable_key(5,fid) /= 0 ) then
       counter = 1
       offset = global_z_indices(START,YBLOCK,myid) - 1
       do k = my_z0(fid)+offset,my_z1(fid)+offset,my_zinc(fid)
        scratch(counter) = s2_bar(k,1)*scalar_scale(2)
        counter = counter + 1
       enddo
       rcode=NF_PUT_VARA_REAL(ncid,s2_barVarID,start1D,count1D,scratch)
       if(rcode.ne.NF_NOERR) print *,myid,'NetCDF ERROR s2_barVarID: ',rcode
      endif      
      
      rcode=NF_CLOSE(ncid)
      if(rcode.ne.NF_NOERR) print *,myid,'ERROR CLOSING NETCDF FILE'
        
999   continue
      call mpi_barrier(comm,ierr)  
      return
     end




subroutine write_netcdf(fid)
 use decomposition_params
 use dimensional_scales
 use dependent_variables
 use independent_variables
 use intermediate_variables, only: div_u,phi,pd
 use methods_params,         only: do_second_scalar
 use mpi_params
 use io_params
 use etc
      
 implicit none
 include 'netcdf.inc'
 

 character(len=80)         :: err_msg
 integer                   :: i,j,k
 integer                   :: Nd_start(4)=0
 integer                   :: Nd_count(4)=0
 integer                   :: ii,jj,kk,offset
 real(kind=4),allocatable  :: scratch(:,:,:)
 real(kind=4)              :: xx_4
     
 integer      :: fid,ncid,rcode  
 integer      :: tVarID,s1VarID,s2VarID   
 integer      :: divustarVarID,phiVarID,pdVarID
 integer      :: uVarID,vVarID,wVarID
 integer      :: do_bar

 if( .NOT. do_write(fid) ) then  !! no data to write
  goto 999
 endif
 
!--------------------------------------------------- 
!  open up the existing, initialized netcdf file
!--------------------------------------------------- 
  rcode=NF_OPEN(fullname(fid),NF_WRITE,ncid)
  if(rcode.ne.NF_NOERR) then
   write(0,*) '... ERROR OPENING NETCDF FILE: write_netcdf ',fid,trim(fullname(fid))
   write(0,*) '... myid, rcode ',myid, rcode
   stop
  endif

!------------------------------------------------------ 
!   extract the variable ids given the variable names
!    only look for ids for variables to be written
!------------------------------------------------------   
  rcode=NF_INQ_VARID(ncid,'time',tVarID)
  if (rcode.ne.NF_NOERR) print *,myid,   &
      'NetCDF ERROR INQ_VARID -> time: ',rcode

  if( variable_key(4,fid) /= 0 ) then
   rcode=NF_INQ_VARID(ncid,'s1',s1VarID)
   if (rcode.ne.NF_NOERR) print *,myid,   &
      'NetCDF ERROR INQ_VARID -> s1: ',rcode
  endif
      
  if( do_second_scalar .and. variable_key(5,fid) /= 0 ) then
   rcode=NF_INQ_VARID(ncid,'s2',s2VarID)
   if (rcode.ne.NF_NOERR) print *,myid,   &
   'NetCDF ERROR INQ_VARID -> s2: ',rcode
  endif
      
  if( variable_key(1,fid) == 1 ) then
   rcode=NF_INQ_VARID(ncid,'u',uVarID)
   if (rcode.ne.NF_NOERR) print *,myid,   &
       'NetCDF ERROR INQ_VARID -> u: ',rcode
  endif
      
  if( variable_key(2,fid) == 1 ) then
      rcode=NF_INQ_VARID(ncid,'v',vVarID)
     if (rcode.ne.NF_NOERR) print *,myid,   &
        'NetCDF ERROR INQ_VARID -> v: ',rcode
  endif
      
  if( variable_key(3,fid) == 1 ) then
     rcode=NF_INQ_VARID(ncid,'w',wVarID)
     if (rcode.ne.NF_NOERR) print *,myid,   &
        'NetCDF ERROR INQ_VARID -> w: ',rcode
  endif
      
  if( variable_key(6,fid) == 1 ) then        
    rcode=NF_INQ_VARID(ncid,'divustar',divustarVarID)
    if (rcode.ne.NF_NOERR) print *,myid,   &
        'NetCDF ERROR INQ_VARID -> divustar: ',rcode
  endif
      
  if( variable_key(7,fid) == 1 ) then        
   rcode=NF_INQ_VARID(ncid,'phi',phiVarID)
   if (rcode.ne.NF_NOERR) print *,myid,   &
       'NetCDF ERROR INQ_VARID -> phi: ',rcode
  endif
      
  if( variable_key(8,fid) == 1 ) then        
   rcode=NF_INQ_VARID(ncid,'pd',pdVarID)
   if (rcode.ne.NF_NOERR) print *,myid,   &
       'NetCDF ERROR INQ_VARID -> pd: ',rcode
  endif
        
!------------------------------------------------------------  
!  set up count and start arrays for this file
!------------------------------------------------------------
  Nd_start(dimid_x(fid))=1
  Nd_count(dimid_x(fid))=count(dimid_x(fid),fid)
       
  Nd_start(dimid_y(fid))=1
  Nd_count(dimid_y(fid))=count(dimid_y(fid),fid)
       
  Nd_start(dimid_z(fid))=1
  Nd_count(dimid_z(fid))=count(dimid_z(fid),fid)
       
  Nd_start(dimid_t(fid))=time_counter(fid)
  Nd_count(dimid_t(fid))=1   !! always 1 time slice written per call
  offset = global_z_indices(START,YBLOCK,myid)-1
  
!------------------------------------------------------------  
!  convert time value to real*4 and write to file
!------------------------------------------------------------  
  xx_4 = t_secs
  rcode=NF_PUT_VARA_REAL(ncid,tVarID,Nd_start(dimid_t(fid)),        &
                         Nd_count(dimid_t(fid)),xx_4)
  if(rcode.ne.NF_NOERR) print *,myid,'NetCDF ERROR tVarID: ',rcode
  if( trim(mode(fid)) == 'append' ) time_counter(fid)=time_counter(fid)+1
  

!------------------------------------------------------------  
!  allocate a contiguous array for selected, 
!  possibly strided output data
!------------------------------------------------------------ 
  allocate( scratch( Nd_count(dimid_x(fid)), &
                     Nd_count(dimid_y(fid)), &
                     Nd_count(dimid_z(fid)) ), stat=rcode )
  if(rcode /= 0 ) stop 'problem allocating scratch in write_netcdf'
  
  
 !------------------------------------------------------------
 ! Scalar 1  in dimensional units.
 !------------------------------------------------------------
  if( variable_key(4,fid) /= 0 ) then
   if( write_s1_bar(fid) ) then
    do_bar=1
   else
    do_bar=0
   endif
   ii=1 ;  
    do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     jj=1;
      do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       kk=1;
        do k=my_z0(fid),my_z1(fid),my_zinc(fid)
          
           scratch(ii,jj,kk) = ( do_bar*s1_bar(k+offset,1)           &
                              + s1(j,i,k))*scalar_scale(1)
          
         kk=kk+1
        enddo
       jj=jj+1
      enddo
     ii=ii+1
    enddo

 !---------------------------------------------------
 !   Each processor writes its local data to 
 !   the appropriate locations in the netcdf file.
 !---------------------------------------------------
   rcode=NF_PUT_VARA_REAL(ncid,s1VarID,Nd_start,Nd_count,scratch)
   if(rcode.ne.NF_NOERR) then
    print *,myid,'NetCDF ERROR: PUT_VARA_REAL scalar (1)',rcode
    err_msg=NF_STRERROR(rcode)
    print*,myid, err_msg
   endif
  endif
      
! Scalar 2  in dimensional units.
  if( do_second_scalar .and. variable_key(5,fid) /= 0 ) then
   if( write_s2_bar(fid) ) then
    do_bar=1
   else
    do_bar=0
   endif
   ii=1 ;  
    do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     jj=1;
      do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       kk=1;
        do k=my_z0(fid),my_z1(fid),my_zinc(fid)
          
           scratch(ii,jj,kk) = ( do_bar*s2_bar(k+offset,1)           &
                              + s2(j,i,k))*scalar_scale(2)
          
         kk=kk+1
        enddo
       jj=jj+1
      enddo
     ii=ii+1
    enddo

 !---------------------------------------------------
 !   Each processor writes its local data to 
 !   the appropriate locations in the netcdf file.
 !---------------------------------------------------
   rcode=NF_PUT_VARA_REAL(ncid,s2VarID,Nd_start,Nd_count,scratch)
   if(rcode.ne.NF_NOERR) then
    print *,myid,'NetCDF ERROR: PUT_VARA_REAL scalar (2)',rcode
    err_msg=NF_STRERROR(rcode)
    print*,myid, err_msg
   endif
  endif
      
! u velocity component in dimensional units.
  if( variable_key(1,fid)==1 ) then
   ii=1 ;  
    do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     jj=1;
      do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       kk=1;
        do k=my_z0(fid),my_z1(fid),my_zinc(fid)
         
          scratch(ii,jj,kk) = u(j,i,k)*velocity_scale
          
          kk=kk+1
        enddo
       jj=jj+1
      enddo
     ii=ii+1
    enddo

!  Each processor writes its local data to 
!  the appropriate locations in the netcdf file.
   rcode=NF_PUT_VARA_REAL(ncid,uVarID,Nd_start,Nd_count,scratch)
   if(rcode.ne.NF_NOERR) print *,myid,   &
      'NetCDF ERROR: PUT_VARA_DOUBLE u',rcode
   endif
      
!      v velocity component in dimensional units.
  if( variable_key(2,fid)==1 ) then
   ii=1 ;  
    do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     jj=1;
      do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       kk=1;
        do k=my_z0(fid),my_z1(fid),my_zinc(fid)
         
          scratch(ii,jj,kk) = v(j,i,k)*velocity_scale
          
          kk=kk+1
        enddo
       jj=jj+1
      enddo
     ii=ii+1
    enddo

!   Each processor writes its local data to 
!   the appropriate locations in the netcdf file.
    rcode=NF_PUT_VARA_REAL(ncid,vVarID,Nd_start,Nd_count,scratch)
    if(rcode.ne.NF_NOERR) print *,myid,   &
       'NetCDF ERROR: PUT_VARA_DOUBLE v',rcode
    endif
      
! w velocity component in dimensional units.
  if( variable_key(3,fid)==1 ) then
   ii=1 ;  
    do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     jj=1;
      do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       kk=1;
        do k=my_z0(fid),my_z1(fid),my_zinc(fid)
         
          scratch(ii,jj,kk) = w(j,i,k)*velocity_scale
                    
          kk=kk+1
        enddo
       jj=jj+1
      enddo
     ii=ii+1
    enddo

!   Each processor writes its local data to 
!   the appropriate locations in the netcdf file.
       
    rcode=NF_PUT_VARA_REAL(ncid,wVarID,Nd_start,Nd_count,scratch)
    if(rcode.ne.NF_NOERR) print *,myid,   &
       'NetCDF ERROR PUT_VARA_DOUBLE: w ',rcode
    endif

! div ustar in dimensional units. ( {1/dt*div(u*)} is in div_u)
  if( variable_key(6,fid)==1 ) then
   ii=1 ;  
    do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     jj=1;
      do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       kk=1;
        do k=my_z0(fid),my_z1(fid),my_zinc(fid)
         
          scratch(ii,jj,kk) = dt*div_u(j,i,k)/(time_scale)
          
         kk=kk+1
        enddo
       jj=jj+1
      enddo
     ii=ii+1
    enddo

!   Each processor writes its local data to 
!   the appropriate locations in the netcdf file.
    rcode=NF_PUT_VARA_REAL(ncid,divustarVarID,Nd_start,Nd_count,scratch)
    if(rcode.ne.NF_NOERR) print *,myid,   &
       'NetCDF ERROR PUT_VARA_DOUBLE: div ustar ',rcode
    endif

! pressure soln in dimensional units.
  if( variable_key(7,fid)==1 ) then
   ii=1 ;  
    do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     jj=1;
      do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       kk=1;
        do k=my_z0(fid),my_z1(fid),my_zinc(fid)
         
          scratch(ii,jj,kk) = phi(j,i,k)*pressure_scale
          
          kk=kk+1
        enddo
       jj=jj+1
      enddo
     ii=ii+1
    enddo

!   Each processor writes its local data to 
!   the appropriate locations in the netcdf file.
    rcode=NF_PUT_VARA_REAL(ncid,phiVarID,Nd_start,Nd_count,scratch)
    if(rcode.ne.NF_NOERR) print *,myid,   &
       'NetCDF ERROR PUT_VARA_DOUBLE: phi ',rcode
    endif
      
! perturbation density in dimensional units.
  if( variable_key(8,fid)==1 ) then
   ii=1 ;  
    do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     jj=1;
      do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       kk=1;
        do k=my_z0(fid),my_z1(fid),my_zinc(fid)
         
          scratch(ii,jj,kk) = pd(j,i,k)*(density_scale)
          
          kk=kk+1
        enddo
       jj=jj+1
      enddo
     ii=ii+1
    enddo

!   Each processor writes its local data to 
!   the appropriate locations in the netcdf file.
    rcode=NF_PUT_VARA_REAL(ncid,pdVarID,Nd_start,Nd_count,scratch)
    if(rcode.ne.NF_NOERR) print *,myid,   &
       'NetCDF ERROR PUT_VARA_REAL: pd ',rcode
    endif

!   Each processor closes the netcdf file.
    rcode=NF_CLOSE(ncid)
    if(rcode.ne.NF_NOERR) print *,myid,   &
      'NetCDF ERROR: closing file in write_netcdf',rcode
    deallocate( scratch )

999   continue
      call mpi_barrier(comm,rcode)  
    return
   end



subroutine process_dimensions(fid,dimstr)
 use io_params
 implicit none
 integer                     :: fid,idim,jdim,ndims
 character(len=3)            :: dimstr
 
 !-------------------------------------------------------------
 ! count the number of spatially varying indices up from 1
 ! if not varying, count down from 4
 ! nspace is the number of spatially varying indices
 ! these values are stored in dimid_x,dimid_y,dimid_z,dimid_t
 !
 ! e.g. all spatial coordinates have > 1 output point
 !      dimid_x,dimid_y,dimid_z,dimid_t ==> 1,2,3,4
 !      nspace = 3
 !
 ! eg.  say y has only 1 value to be output
 !      dimid_x,dimid_z,dimid_t,dimid_y ==> 1,2,3,4
 !      nspace = 2
 !
 ! eg.  say both y and x have only 1 value to be output
 !      dimid_z,dimid_t,dimid_x,dimid_y ==> 1,2,3,4
 !      nspace = 1
 !
 ! eg.  say all of x,y,z have only 1 value to be output
 !      dimid_t,dimid_x,dimid_y,dimid_z ==> 1,2,3,4
 !      nspace = 0
 !-------------------------------------------------------------
  
 idim=1  !! count up for spatially varying dimensions
 jdim=4  !! count down for dimensions w/ indices held constant
 nspace(fid)=0    
  
 if( ilocs(2,fid).ne.ilocs(1,fid) ) then  ! x nonsingleton output coord
  nspace(fid)=nspace(fid)+1
  dimid_x(fid) = idim
  idim = idim + 1
 else
  dimid_x(fid) = jdim
  jdim = jdim - 1
 endif
   
 if( jlocs(2,fid).ne.jlocs(1,fid) ) then  ! y nonsingleton output coord
  nspace(fid)=nspace(fid)+1
  dimid_y(fid) = idim
  idim = idim + 1
 else
  dimid_y(fid) = jdim
  jdim = jdim - 1
 endif
   
 if( klocs(2,fid).ne.klocs(1,fid) ) then  ! z nonsingleton output coord
  nspace(fid)=nspace(fid)+1
  dimid_z(fid) = idim
  idim = idim + 1
 else
  dimid_z(fid) = jdim
  jdim = jdim - 1
 endif
 dimid_t(fid) = idim
 ndims=idim   !! time + counted number of space dims

 !------------------------------------------------------------------
 ! given # of space dimensions, decide which directory to write to
 !------------------------------------------------------------------
 if( ndims-1==0 ) then
  dimstr='TS/'
 elseif( ndims-1==1 ) then
  dimstr='1D/'
 elseif( ndims-1==2 ) then
  dimstr='2D/'
 elseif( ndims-1==3 ) then
  dimstr='3D/'
 endif
end subroutine process_dimensions

subroutine construct_filename(fid,dimstr,topdir,istep,filename)
 use io_params
 use mpi_params,           only: myid,numprocs
 use decomposition_params, only: YBLOCK,proc_row,proc_col
 implicit none
 integer             :: fid
 integer             :: istep
 character(len=3)    :: dimstr
 character(len=80)   :: topdir
 character(len=80)   :: filename
 character(len=6)    :: cnum
 
 character(len=3)    :: c_row_id
 character(len=3)    :: c_col_id
 character(len=7)    :: cid
 

 write(unit = cnum, fmt = 100) istep
 write(unit = c_row_id,  fmt = 101) proc_row(YBLOCK,myid)
 write(unit = c_col_id,  fmt = 101) proc_col(YBLOCK,myid)
 
 cid = c_row_id//'-'//c_col_id
 !write(unit = cid,  fmt = 101) myid


 if( trim(mode(fid)) == 'new' ) then
  if(numprocs .gt. 1) then
   filename = trim(topdir)//dimstr//trim(filename_root(fid))//'_'//cnum//'_'//cid//'.nc'
  else
   filename = trim(topdir)//dimstr//trim(filename_root(fid))//'_'//cnum//'.nc'
  endif
 elseif( trim(mode(fid)) == 'append' ) then
  if(numprocs .gt. 1) then
   filename = trim(topdir)//dimstr//trim(filename_root(fid))//'_'//cid//'.nc'
  else
   filename = trim(topdir)//dimstr//trim(filename_root(fid))//'.nc'
  endif
 endif
 
100 format(I6.6)
101 format(I3.3)
end subroutine construct_filename



subroutine set_local_indices_count(fid)
use io_params
use decomposition_params
use mpi_params,             only: myid
implicit none
integer                        :: fid
integer                        :: offset(2)
integer,external ::  count_vals

 do_write(fid)=.TRUE.   ! start by assuming there is data to write

 !----------------------------------------------------------------
 ! Compute starting and ending indices into local YBLOCK storage
 ! store counts and increments as well
 !----------------------------------------------------------------
 
 !----------------------------------------------------------------
 ! YBLOCK ==> y coord is always local to myid
 !----------------------------------------------------------------
  my_y0(fid)   = jlocs(1,fid)
  my_y1(fid)   = jlocs(2,fid)
  my_yinc(fid) = jlocs(3,fid)
  my_ny(fid)   = count_vals(jlocs(1,fid),jlocs(2,fid),jlocs(3,fid))
  count(dimid_y(fid),fid) = my_ny(fid)
 
 
 !----------------------------------------------------------------
 ! YBLOCK ==> x coord is distributed
 !----------------------------------------------------------------
  offset(1) = global_x_indices(START,YBLOCK,myid)
  offset(2) = global_x_indices(END,YBLOCK,myid)
 
 !-------------------------------------------------------------
 ! if desired data has end index lower than myid's start, or
 ! if desired data has start index higher than myid's end
 ! ===> no data on myid to write for this file set
 !-------------------------------------------------------------
  if( ilocs(1,fid) > offset(2) .or. ilocs(2,fid) < offset(1) ) then
   do_write(fid)=.FALSE.
   return
  endif
 
 !-------------------------------------------------------------
 ! myid has some data to be written, determine local
 ! starting and ending indices and number of points
 !-------------------------------------------------------------
  my_x0(fid)   = maxval( (/ilocs(1,fid), offset(1)/) ) - offset(1) + 1
  my_x1(fid)   = minval( (/ilocs(2,fid), offset(2)/) ) - offset(1) + 1 
  my_xinc(fid) = ilocs(3,fid)
  my_nx(fid)   = count_vals(my_x0(fid),my_x1(fid),my_xinc(fid))
  count(dimid_x(fid),fid) = my_nx(fid)
  
  
 !----------------------------------------------------------------
 ! YBLOCK ==> z coord is distributed
 !----------------------------------------------------------------
  offset(1) = global_z_indices(START,YBLOCK,myid)
  offset(2) = global_z_indices(END,YBLOCK,myid)
 
 !-------------------------------------------------------------
 ! if desired data has end index lower than myid's start, or
 ! if desired data has start index higher than myid's end
 ! ===> no data on myid to write for this file set
 !-------------------------------------------------------------
  if( klocs(1,fid) > offset(2) .or. klocs(2,fid) < offset(1) ) then
   do_write(fid)=.FALSE.
   return
  endif
 
 !-------------------------------------------------------------
 ! myid has some data to be written, determine local
 ! starting and ending indices and number of points
 !-------------------------------------------------------------
  my_z0(fid)   = maxval( (/klocs(1,fid), offset(1)/) ) - offset(1) + 1
  my_z1(fid)   = minval( (/klocs(2,fid), offset(2)/) ) - offset(1) + 1 
  my_zinc(fid) = klocs(3,fid)
  my_nz(fid)   = count_vals(my_z0(fid),my_z1(fid),my_zinc(fid)) 
  count(dimid_z(fid),fid) = my_nz(fid)
   
  return 
 end subroutine set_local_indices_count
 
 integer function count_vals(i0,i1,inc) result(numvals)
  integer, intent(in)  :: i0,i1,inc
  if( i1<i0 ) then
   numvals=0
   return
  endif
  numvals = floor( float(i1-i0)/inc ) + 1
end function count_vals





subroutine write_particle_data
  use particles
  use independent_variables, only: t_secs
  use dimensional_scales,    only: length_scale
  use mpi_params,            only: myid,comm,numprocs
      
  implicit none
  include 'netcdf.inc'
  
  logical, SAVE                 :: first_entry=.TRUE.
  integer,save                  :: iid  !!dimension id for particles
  integer,save                  :: jid  !!dimension id for time
  integer,save                  :: idimensionSIZE(1) 
  integer,save                  :: jdimensionSIZE(1) 
  integer,save                  :: idimensionID !! for dimension variable "float id"
  integer,save                  :: tVarID !! for dimension variable "time"
  integer,save                  :: rcode, ierr
  integer,save                  :: ncid
      
  !!Variables: x(particle id, time),
  !!           y(particle id, time),
  !!           z(particle id, time)
  integer,save                  :: All2DVarSIZE(2),start2D(2),count2D(2)
  integer,save                  :: xVarID,yVarID,zVarID   !! 2D Variable ids
  integer,save                  :: istart,jstart,icount,jcount,i
  integer,allocatable,save      :: idata(:)
  real(kind=4),allocatable,save :: xdata(:)
  real(kind=4),save             :: xx_4
  character(len=80),save        :: topdir
  character(len=80),save        :: filename
  character(len=80),save        :: filename_root
  character(len=4),save         :: cid
  logical,save                  :: COMPAS_SCRATCH = .FALSE.
    
    
  if( first_entry ) then
    
    j_particle_time = 1
    filename_root='particles'
    write(unit = cid,  fmt = 101) myid
101 format(I4.4)   
    if( COMPAS_SCRATCH ) THEN
     topdir='/scratch/output/lagrangian/'     
    else
     topdir='output/lagrangian/'
    endif
   
    if(numprocs .gt. 1) then
     filename = trim(topdir)//trim(filename_root)//'_'//cid//'.nc'
    else
     filename = trim(topdir)//trim(filename_root)//'.nc'
    endif
     
  
   
   icount = my_last - my_first + 1
   allocate( idata(icount), stat=rcode )
   if( rcode /= 0 ) &
    stop 'problem allocating idata in write_particle_data'
   allocate( xdata(icount), stat=rcode  )
   if( rcode /= 0 ) &
    stop 'problem allocating idata in write_particle_data'
 
   !!write(0,*) 'creating particle file ',t_secs,init_particle_file
   rcode = NF_CREATE(filename, NF_NOCLOBBER, ncid) 
   if(rcode.ne.NF_NOERR) then
    write(0,*) 'ERROR CREATING NETCDF FILE FOR PARTICLES ',filename
    stop
   endif

   !!Dimensions Definition
   icount = my_last - my_first + 1
   rcode = NF_DEF_DIM(ncid, 'idimension', icount, iid)
   if(rcode.ne.NF_NOERR) write(0,*) 'init for particles NetCDF ERROR: ',rcode

   rcode = NF_DEF_DIM(ncid, 'jdimension', NF_UNLIMITED, jid)
   if(rcode.ne.NF_NOERR) print *,'init for particles NetCDF ERROR: ',rcode

   !!Initialize size arrays for dimensions
   idimensionSIZE(1) = iid
   jdimensionSIZE(1) = jid

   !!Define dimension variables
   rcode = NF_DEF_VAR(ncid,'pid',NF_INT,1,idimensionSIZE,idimensionID)
   if(rcode.ne.NF_NOERR) print *,'NetCDF ERROR: ',rcode
   rcode=NF_PUT_ATT_TEXT(ncid,idimensionID,'long_name',11,'particle id')

   rcode = NF_DEF_VAR(ncid,'t_secs',NF_FLOAT,1,jdimensionSIZE,tVarID)
   if(rcode.ne.NF_NOERR) print *,'NetCDF ERROR: ',rcode
   rcode=NF_PUT_ATT_TEXT(ncid,tVarID,'long_name',4,'time')
   rcode=NF_PUT_ATT_TEXT(ncid,tVarID,'units',7,'seconds')

   !! Define 2D position variables
   !!Initialize size arrays for id/time varying variables
   All2DVarSIZE(1) = iid
   All2DVarSIZE(2) = jid


   !! x particle position = x
   rcode = NF_DEF_VAR(ncid,'x',NF_FLOAT,2,All2DVarSIZE,xVarID)
   if(rcode.ne.NF_NOERR) print *,'NetCDF ERROR: ',rcode
   rcode=NF_PUT_ATT_TEXT(ncid,xVarID,'long_name',19,'x particle position')                      
   rcode=NF_PUT_ATT_TEXT(ncid,xVarID,'units',1,'m')
  
   !! y particle position = y
   rcode = NF_DEF_VAR(ncid,'y',NF_FLOAT,2,All2DVarSIZE,yVarID)
   if(rcode.ne.NF_NOERR) print *,'NetCDF ERROR: ',rcode
   rcode=NF_PUT_ATT_TEXT(ncid,yVarID,'long_name',19,'y particle position')                      
   rcode=NF_PUT_ATT_TEXT(ncid,yVarID,'units',1,'m')
  
   !! z particle position = z
   rcode = NF_DEF_VAR(ncid,'z',NF_FLOAT,2,All2DVarSIZE,zVarID)
   if(rcode.ne.NF_NOERR) print *,'NetCDF ERROR: ',rcode
   rcode=NF_PUT_ATT_TEXT(ncid,zVarID,'long_name',19,'z particle position')                      
   rcode=NF_PUT_ATT_TEXT(ncid,zVarID,'units',1,'m')
     
   !!End define mode
   rcode=NF_ENDDEF(ncid)
   if(rcode.ne.NF_NOERR) print *,'ERROR  LEAVING DEFINE MODE'
  

   !!Write out dimension variable: particle id
   do i = 1,icount
    idata(i) = my_first + i - 1   !!--> integer
   enddo

   !! Put idimension data, i.e. particle id, into netcdf file
   istart = 1
   rcode=NF_PUT_VARA_INT(ncid,idimensionID,istart,icount,idata)
   if(rcode.ne.NF_NOERR) print *,'NetCDF ERROR storing particle ids: ',rcode

   !!Close netcdf file
   rcode=NF_CLOSE(ncid)
   if(rcode.ne.NF_NOERR) print *,'ERROR CLOSING particle NETCDF FILE'
  
   !!Initialize size arrays
   start2D(1) = 1
   start2D(2) = j_particle_time
   count2D(1) = icount
   count2D(2) = 1
  
  first_entry=.FALSE.
endif  ! end file initialization 



!! Open the previously created netcdf data file
  rcode=NF_OPEN(filename,NF_WRITE,ncid)
  if(rcode.ne.NF_NOERR) then 
   print *,'ERROR OPENING NETCDF FILE: write_netcdf ',filename,ncid
   stop
  endif
   
  start2D(2) = j_particle_time
  count2D(2) = 1
  
!!Write out time in seconds
  jstart = j_particle_time
  jcount = 1
  rcode=NF_INQ_VARID(ncid,'t_secs',tVarID)
  if(rcode.ne.NF_NOERR) print *,'NetCDF ERROR: NF_INQ_VARID ',rcode
  
  xx_4=t_secs  !! convert to float
  rcode=NF_PUT_VARA_REAL(ncid,tVarID,jstart,jcount,xx_4)
  if(rcode.ne.NF_NOERR) print *,'NetCDF ERROR: NF_PUT_VARA_REAL',rcode


  !!Store x float position variable x
  do  i=1,count2D(1)
   xdata(i) = positions(my_first+i-1,1)*length_scale
  enddo
  rcode=NF_INQ_VARID(ncid,'x',xVarID)
  rcode=NF_PUT_VARA_REAL(ncid,xVarID,start2D,count2D,xdata)
  if(rcode.ne.NF_NOERR) print *,'NetCDF ERROR: ',rcode

  !!Store y float position variable y
  do  i=1,count2D(1)
   xdata(i) = positions(my_first+i-1,2)*length_scale
  enddo
  rcode=NF_INQ_VARID(ncid,'y',yVarID)
  rcode=NF_PUT_VARA_REAL(ncid,yVarID,start2D,count2D,xdata)
  if(rcode.ne.NF_NOERR) print *,'NetCDF ERROR: ',rcode

  !!Store z float position variable z
  do  i=1,count2D(1)
   xdata(i) = positions(my_first+i-1,3)*length_scale
  enddo
  rcode=NF_INQ_VARID(ncid,'z',zVarID)
  rcode=NF_PUT_VARA_REAL(ncid,zVarID,start2D,count2D,xdata)
  if(rcode.ne.NF_NOERR) print *,'NetCDF ERROR: ',rcode

  !! Close netcdf file
  rcode=NF_CLOSE(ncid)
  if(rcode.ne.0) print *,'ERROR CLOSING NETCDF FILE'
  
  !!increment integer time index
  j_particle_time = j_particle_time + 1 

 call mpi_barrier(comm,ierr)  !! sync up here
 return
end subroutine write_particle_data


