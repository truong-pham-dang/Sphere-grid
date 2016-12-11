!  Sphere_grid_fortran.f90 
!
!  FUNCTIONS:
!  Sphere_grid_fortran - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Sphere_grid_fortran
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

!*******************************************************************************
! Grid generation code for a sphere.
!
! This is Version 4 (April 11, 2016).
!
! This code generates a 3D grid over a sphere centered at the origin.
!
! Note: Use ifort to compile this code. Other compilers may not work
!       for the binary-file writing (for .ugrid). Not sure how to fix it...
!
! ------------------------------------------------------------------------------
!
!  Input:
!
!                    rmax = thickness of a prismatic layer in mixed grids.
!  target_reynolds_number = Reynolds number to determine the mesh spacing at wall
!           target_y_plus = y-plus to determine the mesh spacing at wall
!                           or a user-specified spacing if the input Re < 0.
!
!              igrid_type = Element type. 1=Prism, 2=Tet, 3=Mixed(tet/prism)
!         b8_ugrid_format = T: unformatted(binary file), F: formatted 
!                distance = Distance to the outer boundary from the origin.
!                   nr_gs = Number of elements from the apex to the shoulder
! Output:
!
!   Two output files by default:
!
!     "sphere_'element_type'.1.mapbc"            !Boundary condition file for FUN3D
!     "sphere_'element_type'.1.ugrid/.b8.ugrid"  !Unstructured grid
!     "sphere_'element_type'_boundary_tec.1.dat" !Boundary grid data for Tecplot
!
!  --------
!   Below are optional (can be activated inside the program):
!
!   (1)Tecplot files for viewing the 3D grid:
!
!      "sphere_'element_type'_volume_tec.1.dat"
!
!      It can be activated by setting generate_tec_file_v = .true.
!
!   (2)Line information
!
!     "sphere_'element_type'lines_fmt"   !Lines for implicit solvers or line-agglomeration.
!
!     It can be activated by setting generate_line_file = .true.
!
!   (3)File containing k1-k2-k3-k4 structured indices.
!
!     "sphere_'element_type'k"
!
!     It can be activated by setting generate_k_file = .true.
!     Set line_extent also to select the type of lines:
!      1 = Lines within a thin boundary layer region
!      2 = Lines go all the way to the outer boundary
!
!
!        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!
! the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!
! Version 4
! 04-11-2016: Comments updated.
!
! This F90 code is written and made available for an educational purpose.
! This file may be updated in future.
!
! Katate Masatsuka, http://www.cfdbooks.com
!
!*******************************************************************************
 program sphere_grid

 implicit none

! Parameters
  integer , parameter ::    dp = selected_real_kind(P=15)
  real(dp), parameter ::  zero = 0.0_dp
  real(dp), parameter ::   one = 1.0_dp
  real(dp), parameter ::   two = 2.0_dp
  real(dp), parameter :: three = 3.0_dp
  real(dp), parameter ::   six = 6.0_dp
  real(dp), parameter ::  half = 0.5_dp
  real(dp), parameter ::    pi = 3.14159265358979323846_dp

! Custom data types
  type tria_data
   integer, dimension(3) :: v    !Vertices (nodes) of the triangle
   integer               :: type !Type of triangle (upward, downward, cylinder)
  end type tria_data

  type node_data_xy
    integer :: gnode     !Global node number
   real(dp) :: x, y, z   !Coordinates in xy-plane.
  end type node_data_xy

  type node_data
   real(dp) :: x,y,z       !Coordinates in xyz space
   real(dp) :: nx,ny,nz    !Unit vector normal to the surface
   real(dp) :: nx2,ny2,nz2 !Unit vector along the second direction
   logical  :: circum_node
  end type node_data

! Input variables
   integer :: igrid_type !Type of grid: 1=Prismatic, 2=Tet, 3=Mixed
   integer :: nr_gs      !Number of division of the generating sector

! Output File names

! (1) Tecplot files for viewing intermediate/final grids for debugging.
  character(80) :: filename_gs        = "debug_generating_sector.dat"
  character(80) :: filename_disk      = "debug_disk.dat"
  character(80) :: filename_surface   = "debug_hemisphere_surface.dat"
  character(80) :: filename_surface_s = "debug_sphere_surface.dat"
  character(80) :: filename_tecplot_b
  character(80) :: filename_tecplot_v

! (2) These are the files used for CFD computations
  character(80) :: filename_mapbc
  character(80) :: filename_lines
  character(80) :: filename_ugrid
  character(80) :: filename_k

! Local variables
   integer :: os           !IO constant
  real(dp) :: distance     !Distance to the outer boundary from the body
   integer :: i,j,k,inode

   integer :: ntrias_gs    !Number of triangles in the generating sector
   integer :: nnodes_gs    !Number of nodes in the generating sector
  real(dp) :: Rd           !Radius of the sphere
  real(dp) :: r_gs         !Radius of the generating sector
  real(dp) :: dr_gs        !Spacing along the side of the generating sector
  real(dp) :: r2, dtheta, theta
  type(node_data_xy), dimension(:),     pointer :: nodes1, nodes2, node_gs
  integer, dimension(:),     pointer :: k1_gs, k2_gs, k3_gs
  type(tria_data)   , dimension(:), allocatable :: tria_gs

  integer :: nnodes_disk, ntrias_disk
  type(node_data_xy), dimension(:),     pointer :: node_disk
  integer, dimension(:),     pointer :: k1_disk
  integer, dimension(:),     pointer :: k2_disk
  integer, dimension(:),     pointer :: k3_disk
  type(tria_data)   , dimension(:), allocatable :: tria_disk

  real(dp) :: xp, yp, zp

  integer  :: nnodes_hs
  real(dp) :: s
  type(node_data), dimension(:),     pointer :: node_body_hs
  type(tria_data), dimension(:), allocatable :: tria_hs
  integer, dimension(:),     pointer :: k1_body_hs, k2_body_hs, k3_body_hs, k4_body_hs
  integer :: ntrias_hs

  integer  :: nnodes
  type(node_data), dimension(:),     pointer :: node_body
  type(tria_data), dimension(:), allocatable :: tria
  integer, dimension(:),     pointer :: k1_body, k2_body, k3_body, k4_body
  integer :: ntrias

  integer, dimension(:), allocatable :: node_map
  integer :: node1, node2, node3

  integer :: nnodes_circum
  integer, dimension(:)  , allocatable :: nodes_circum
  type(node_data), dimension(:), pointer :: node
  integer, dimension(:)  , allocatable :: k1, k2, k3, k4

  integer  :: nr, nm
  real(dp) :: rmax, drN, dr1, dr_outer, gr
  real(dp) :: sf, drnmp1 !Stretching factor in the outer region
  real(dp), allocatable, dimension(:) :: rs
  real(dp), allocatable, dimension(:) :: vspacing
  real(dp), allocatable, dimension(:) :: vspacing_nm

  integer :: nnodes_body
  integer :: node4, node5, node6
  integer :: ntrias_b, nquads_b, ntet, nprs
  integer,  dimension(:,:), allocatable :: tri, quad, tet, prs
  integer,  dimension(:)  , allocatable :: node_above

  real(dp) :: xc, yc, zc
  real(dp) :: dirx, diry, dirz

! Remove 6 nodes (Optional)
 
  integer, dimension(6)   :: nodes_removed, nodes_removed_nghbr
  integer, dimension(6,2) :: nodes_changed
  integer :: rnode, newnode
  logical :: deleted
  character(80) :: remove_corner_nodes

! For smoothing applied on the disk triangulation
  integer  :: i_smoothing, max_smoothing, v(3)
  real(dp) :: smoothing_factor, dxy_norm, dxyz_norm
  real(dp), dimension(:,:), allocatable :: dxy
  real(dp), dimension(:,:), allocatable :: dxyz
  integer  :: smoothing_type

  integer :: line_extent, mglevels

  logical :: b8_ugrid_format     = .false.
  logical :: generate_tec_file_v = .true.
  logical :: generate_tec_file_b = .false.
  logical :: generate_k_file     = .false.
  logical :: generate_line_file  = .false.

  real(dp) :: target_reynolds_number, target_yplus, cf
  logical  :: debug_mode
  integer  :: sk12p, sk23p, sk31p
  integer  :: sk12m, sk23m, sk31m
  integer  :: n_sections

  real(dp) :: ds, ds_outer, rs_nmm1, rs_nmp1, rs_nm, rs_nr, rs_nrm1
  real(dp) :: dr_spacing, err_dr_spacing
  integer  :: i_nr

      nprs = 0
      ntet = 0
  ntrias_b = 0
  nquads_b = 0

   debug_mode = .false.

    n_sections = 6 ! Full geometry, always.

!*******************************************************************************
! Sphere Geometry:
!
!
!        z
!        |                      y
!        |                      / 
!        |           Sphere    /     
!        |                * * /
!        |             *       *
!        |            *     /   *
!        |-----------------o---------------------------------> x
!                     *  Origin * 
!                      *       *
!                         * * 
!
!
!
!  The sphere is centered at the origin with the radius of 0.5.
!
!  Outer boundary is also a spherical surface centered at the origin with the radius
!  s[ecify by the input: 'distance'.
!
!*******************************************************************************

   Rd =  0.5_dp ! Sphere radius

!*******************************************************************************
! 0. Select the type of grid to generate.
!*******************************************************************************
!if (1==0) then

!-------------------------------------------------------------------
   write(*,*) 
   write(*,*) " Assumed boundary layer thickness: rmax = ? (E.g.: Try 0.1)"
   write(*,*) 
   write(*,*) "  (Note: Distance from the wall within which you assume"
   write(*,*) "         the boundary layer exists. This is the region "
   write(*,*) "          where the prismatic layer is constructed in"
   write(*,*) "          case of mixed grids. Nodes will be placed"
   write(*,*) "         normal to the wall by geometric stretching in this region."
    read(*,*) rmax

!-------------------------------------------------------------------
   write(*,*) 
   write(*,*) " To determine the grid spacing above the wall..."
   write(*,*) "  1. Target Reynolds number (to specify y+) = ? (E.g.: Try 10000)"
   write(*,*) 
   write(*,*) "  (Note: Input a negative value to specify the grid spacing in the grid unit. "
   write(*,*) "         E.g., type in -10, then input target_y_plus manually next.         ) "
    read(*,*) target_reynolds_number

   if (target_reynolds_number > zero) then
    write(*,*) 
    write(*,*) "  2. Target y_plus value = ? (E.g.: Try 0.01)"
   else
    write(*,*) 
    write(*,*) "  2. First vertical spacing off the wall = ?"
    write(*,*) "  (Note: The diameter of the sphere is 1.0.)"
   endif

     read(*,*) target_yplus

!-------------------------------------------------------------------
   write(*,*)
   write(*,*) 
   write(*,*) "Grid type = 1: Prismatic grid"
   write(*,*) "          = 2: Tetrahedral grid"
   write(*,*) "          = 3: Mixed grid: Tet/prism"
 1 write(*,*) "grid_type = ? (E.g.: Try 2)"

   read(*,*) igrid_type
   if (igrid_type/=1 .and. igrid_type/=2 .and. igrid_type/=3) then
    write(*,*) " >>> Invalid input value."
    go to 1
   endif

!-------------------------------------------------------------------
   write(*,*) "Use b8.ugrid format for grid (T/F)= ? (E.g.: Try F)"
   write(*,*) " F -> .ugrid    (formatted                   )"
   write(*,*) " T -> .b8.ugrid (unformatted stream bigEndian)"
   read(*,*) b8_ugrid_format

   if ( big_endian_io(9999) ) then
     write(*,*) 'The system is big Endian'
     write(*,*) ' Ensure big Endian -> setenv F_UFMTENDIAN big'
   else
     write(*,*) 'The system is little Endian'
   endif

!-------------------------------------------------------------------
  !(1) Prismatic grid
  if     ( igrid_type ==1  ) then

   filename_mapbc     = "sphere_prism.1.mapbc"
   filename_lines     = "sphere_prism.1.lines_fmt"
   if ( b8_ugrid_format ) then
     if ( big_endian_io(9999) ) then
      filename_ugrid     = "sphere_prism.1.b8.ugrid"
     else
      filename_ugrid     = "sphere_prism.1.lb8.ugrid"
     end if
   else
   filename_ugrid     = "sphere_prism.1.ugrid"
   endif
   filename_k         = "sphere_prism.1.k"
   filename_tecplot_b = "sphere_prism_boundary_tec.1.dat"
   filename_tecplot_v = "sphere_prism_volume_tec.1.dat"

!-------------------------------------------------------------------
  !(2) Tetrahedral grid
  elseif ( igrid_type ==2  ) then

   filename_mapbc     = "sphere_tetra.1.mapbc"
   filename_lines     = "sphere_tetra.1.lines_fmt"
   if ( b8_ugrid_format ) then
     if ( big_endian_io(9999) ) then
      filename_ugrid     = "sphere_tetra.1.b8.ugrid"
     else
      filename_ugrid     = "sphere_tetra.1.lb8.ugrid"
     end if
   else
   filename_ugrid     = "sphere_tetra.1.ugrid"
   endif
   filename_k         = "sphere_tetra.1.k"
   filename_tecplot_b = "sphere_tetra_boundary_tec.1.dat"
   filename_tecplot_v = "sphere_tetra_volume_tec.1.dat"

!-------------------------------------------------------------------
  !(3) Mixed grid: Prism near the surface and Tetra otherwise.
  elseif ( igrid_type ==3  ) then

   filename_mapbc     = "sphere_mixed.1.mapbc"
   filename_lines     = "sphere_mixed.1.lines_fmt"
   if ( b8_ugrid_format ) then
     if ( big_endian_io(9999) ) then
      filename_ugrid     = "sphere_mixed.1.b8.ugrid"
     else
      filename_ugrid     = "sphere_mixed.1.lb8.ugrid"
     end if
   else
   filename_ugrid     = "sphere_mixed.1.ugrid"
   endif
   filename_k         = "sphere_mixed.1.k"
   filename_tecplot_b = "sphere_mixed_boundary_tec.1.dat"
   filename_tecplot_v = "sphere_mixed_volume_tec.1.dat"

  endif

!-------------------------------------------------------------------
!  smoothing_type = Smooth the surface triangulation.
!     0 = no smoothing; 1 smoothing; 2 constrained smoothing
!   write(*,*) 
!   write(*,*) " Smoothing the sphere surface triangulation"
!   write(*,*) "   0 : No smoothing"
!   write(*,*) "   1 : Smoothing"
!   write(*,*) "   2 : Constrained smoothing"
!   read(*,*) smoothing_type

   smoothing_type = 1

!-------------------------------------------------------------------
!     remove_corner_nodes = Remove corner nodes for smooth triangulation.
!                           "yes" = remove;  "no" = don't remove
!   write(*,*) 
!   write(*,*) "--- Node removal recommended if smoothing applied..."
!   write(*,*) " Remove_corner_nodes? (Input yes or no)"
!   write(*,*) " NOTE: This should be 'no' for k1-k2-k3-k4 structured grid"
!   read(*,*) remove_corner_nodes

   remove_corner_nodes = "no"

!-------------------------------------------------------------------
!  Distance to the outer boundary
   write(*,*) 
   write(*,*) "--- Location of outer boundary..."
   write(*,*) " Distance to the outer boundary from the body = ? (E.g., Try 100.0)"
   read(*,*) distance

!-------------------------------------------------------------------
! (Optional) Implicit line type

!   write(*,*) 
!   write(*,*) "--- Implicit lines..."
!   write(*,*) "   1 = Lines within a thin boundary layer region"
!   write(*,*) "   2 = Lines go all the way to the outer boundary."
!   write(*,*) " line_extent = ?"
!   read(*,*) line_extent

    line_extent = 2

!-------------------------------------------------------------------
! (Optional) Optional files

!   write(*,*) "generate_tec_file for volume (T/F) = ?"
!   read(*,*)   generate_tec_file_v
    generate_tec_file_v = .false.

!   write(*,*) "generate_tec_file for boundary (T/F) = ?"
!   read(*,*)   generate_tec_file_b
    generate_tec_file_b = .true.

!   write(*,*) "generate_k_file (T/F) = ?"
!   read(*,*)   generate_k_file
    generate_k_file     = .false.

!   write(*,*) "generate_line_file (T/F) = ?"
!   read(*,*)   generate_line_file
    generate_line_file= .false.

   write(*,*) "***********************************************************"
   write(*,*) " 1. Triangulate the generating sector, the building block"
   write(*,*) "***********************************************************"
   write(*,*) "Division of the generating sector = ? (E.g., Try 6)"

   read(*,*) nr_gs

!*******************************************************************************
! 1. Systematic triangulation of the generating sector (isotropic).
!    Resolution of this sector will determine other dimensions.
!
! This is a sector with the central angle 60 degrees.
! It is called here the generating sector.
! It is located in the xy-plane.
! We first triangulate this, and use it to build a triangulation of the
! whole circle (the disk).
!
!       ____________
!       \/\/\/\/\/\/    This is an example corresponding to nr_gs = 6.
!        \/\/\/\/\/
!         \/\/\/\/        y ^
!          \/\/\/           |
!           \/\/            |
!            \/             ----> x
!
! NOTE: The number of triangles = 1 + 3 + 5 + ... + 2*nr_gs-1 = nr_gs*nr_gs
!       The number of nodes  = 1 + 2 + 3 + ... + nr_gs+1
!
! NOTE: Important to distinguish two types of triangles.
!       They will define the way the prism is subdivided into 3 tets.
!       The subdivision must be done carefully to match the triangular faces
!       between two adjacent tetrahedra.
!
!*******************************************************************************

   k = 0
   do
    k = k + 1
    if (mod(nr_gs,2**(k-1))==0) then
      mglevels = k
    else
      write(*,*)
      write(*,*) " Maximum multigrid level is ", mglevels
      write(*,*)
      exit
    endif
   end do

   r_gs = half*pi*Rd          ! Radius of the isotropic triangle
   write(*,*) " Radius of the isotropic triangle, r_gs = ", r_gs

  dr_gs = r_gs/real(nr_gs,dp) ! Uniform refinement
   write(*,*) " dr_gs = ", dr_gs

  nnodes_gs = (nr_gs+1)*((nr_gs+1)+1)/2
   write(*,*) " nnodes_gs = ", nnodes_gs
  ntrias_gs = nr_gs**2
   write(*,*) " ntrias_gs = ", ntrias_gs
  allocate(node_gs(nnodes_gs))
  allocate(tria_gs(ntrias_gs))
  allocate(k1_gs(nnodes_gs))
  allocate(k2_gs(nnodes_gs))
  allocate(k3_gs(nnodes_gs))

  nnodes_gs = 0
  ntrias_gs = 0

  triangulate : do i = 1, nr_gs

!   r1 = dr_gs * real(i-1,dp) ! Current arc
    r2 = dr_gs * real(i,dp)   !    Next arc

! Nodes on the current arc (r = r1)
   call my_alloc_ndxy_ptr(nodes1,i)

   if (i ==  1) then

      nodes1(1)%x = zero
      nodes1(1)%y = zero
      nodes1(1)%gnode = 1

        nnodes_gs = 1
     node_gs(1)%x = zero
     node_gs(1)%y = zero

     k1_gs(1) = 0
     k2_gs(1) = 0
     k3_gs(1) = - ( k1_gs(1) + k2_gs(1) )

   else

    do k = 1, i
     nodes1(k)%x     = nodes2(k)%x
     nodes1(k)%y     = nodes2(k)%y
     nodes1(k)%gnode = nodes2(k)%gnode
    end do

   endif

! Nodes on the next arc (r = r2): New nodes

    call my_alloc_ndxy_ptr(nodes2,i+1)
    dtheta = (pi/three) / real(i,dp)
    do k = 1, i+1
     theta = dtheta * real(k-1,dp)
     nodes2(k)%x = r2 * cos(theta)
     nodes2(k)%y = r2 * sin(theta)

     nnodes_gs = nnodes_gs + 1
     nodes2(k)%gnode = nnodes_gs
     node_gs(nnodes_gs)%x = nodes2(k)%x
     node_gs(nnodes_gs)%y = nodes2(k)%y


     k1_gs(nnodes_gs) = i + (1 - k)
     k2_gs(nnodes_gs) = k - 1
     k3_gs(nnodes_gs) = - ( k1_gs(nnodes_gs) + k2_gs(nnodes_gs) )


!    Keep the record of corner nodes and their neighbors for node removal
     if (trim(remove_corner_nodes) == "yes") then
      if (i==nr_gs-1 .and. k==1  ) nodes_removed_nghbr(1) = nnodes_gs
      if (i==nr_gs   .and. k==1  ) nodes_removed(1)       = nnodes_gs
      if (i==nr_gs-1 .and. k==i+1) nodes_removed_nghbr(2) = nnodes_gs
      if (i==nr_gs   .and. k==i+1) nodes_removed(2)       = nnodes_gs
     endif

    end do

! Triangulate the region between nodes1 and nodes2.
! Right to left

! NOTE: Nodes are ordered clockwise at this point.
!       It will be switched to counter-clockwise later.

! Type 1 triangle
!
! nodes2(k+1)   nodes2(k)
!      1             2
!       o-----------o
!        \         /  
!         \       /
!          \     /
!           \   /
!            \ /
!             o
!             3
!         nodes1(k)

   downward_tria : do k = 1, i
    ntrias_gs = ntrias_gs + 1
    tria_gs(ntrias_gs)%v(1) = nodes2(k+1)%gnode
    tria_gs(ntrias_gs)%v(2) = nodes2(k  )%gnode
    tria_gs(ntrias_gs)%v(3) = nodes1(k  )%gnode
    tria_gs(ntrias_gs)%type = 1
   end do downward_tria


! Type 2 triangle
!
!         nodes2(k+1)
!             3
!             o
!            / \
!           /   \
!          /     \ 
!         /       \
!        /         \
!       o-----------o
!      2             1
! nodes1(k+1)   nodes1(k)

   if (i > 1) then
    upward_tria : do k = 1, i-1
     ntrias_gs = ntrias_gs + 1
     tria_gs(ntrias_gs)%v(1) = nodes1(k  )%gnode
     tria_gs(ntrias_gs)%v(2) = nodes1(k+1)%gnode
     tria_gs(ntrias_gs)%v(3) = nodes2(k+1)%gnode
     tria_gs(ntrias_gs)%type = 2

     if (trim(remove_corner_nodes) == "yes") then
      if (i==nr_gs .and. k==i-1) then
       tria_gs(ntrias_gs)%type = 20
      endif
     endif

    end do upward_tria
   endif

  end do triangulate

! So, now I have a triangulation defined by tria_gs and node_gs.
! Number of triangles = ntrias_gs
! Number of nodes     = nnodes_gs

!*******************************************************************************
! Write a Tecplot file for the generating sector
!******************************************************************************
 debug_mode_01 : if (debug_mode) then
 open(unit=1, file=filename_gs, status="unknown", iostat=os)

  write(1,*) 'TITLE = "GRID"'
  write(1,*) 'VARIABLES = "x","y","z","k1+k2"'
  write(1,*) 'ZONE  N=', nnodes_gs,',E=', ntrias_gs,' , ET=triangle, F=FEPOINT'

! Nodes
  do i = 1, nnodes_gs
    write(1,'(3ES20.10,i10)') node_gs(i)%x, node_gs(i)%y, 0.0, k1_gs(i)+k2_gs(i)
  end do

! Triangles
  do i = 1, ntrias_gs
   write(1,'(3I10)') tria_gs(i)%v(1), tria_gs(i)%v(2), tria_gs(i)%v(3)
  end do

 close(1)

 write(*,*) "Tecplot file has been written: ", filename_gs
 endif debug_mode_01
!*******************************************************************************

!*******************************************************************************
! 2. Generate a triangulated disk.
!
! Now, rotate and copy the sector triangulation (as indicated by 'TS' below) onto
! 5 places (indicated by 1,2,3,4,5 below) to form a triangulation of a whole disk
! (6 patches in total):
!
!           .  . 
!        \        /
!     .   \   1  /   .
!    .     \    /     .
!   .   2   \  /  TS   .
!   _________\/_________
!   .        /\        .
!   .   3   /  \   5   .
!    .     /    \     .
!     .   /      \   .
!        /   4    \ .
!           .  .
!
! Generate new data: tria_disk and node_disk
! Number of triangles = ntrias_disk
! Number of nodes     = nnodes_disk
!*******************************************************************************

  write(*,*) "***********************************************************"
  write(*,*) " 2. Generate a disk by copying the generating sector"
  write(*,*) "***********************************************************"

  nnodes_disk = nnodes_gs * n_sections !<- More than enough (by the # of overlapping ndoes)
  ntrias_disk = ntrias_gs * n_sections
  allocate(node_disk(nnodes_disk))
  allocate(tria_disk(ntrias_disk))
  allocate(k1_disk(nnodes_disk))
  allocate(k2_disk(nnodes_disk))
  allocate(k3_disk(nnodes_disk))

  nnodes_disk = 0
  ntrias_disk = 0

! Copy the data from the generating-sector data to the disk data.

!  Copy the node data
!  This is the first one: theta = 0 to 60 degrees

   do i = 1, nnodes_gs
    node_disk(i)%x = node_gs(i)%x
    node_disk(i)%y = node_gs(i)%y

        k1_disk(i) = k1_gs(i)
        k2_disk(i) = k2_gs(i)
        k3_disk(i) = k3_gs(i)

   end do

   nnodes_disk = nnodes_gs

!  Copy the triangle data
!  This is the first one: theta = 0 to 60 degrees

   do i = 1, ntrias_gs
    tria_disk(i)%v    = tria_gs(i)%v
    tria_disk(i)%type = tria_gs(i)%type
   end do

   ntrias_disk = ntrias_gs

! Now generate other parts of the disk: i=1,2,3,4,5
! 1. theta =  60 to 120 degrees
! 2. theta = 120 to 180 degrees
! 3. theta = 180 to 240 degrees
! 4. theta = 240 to 300 degrees
! 5. theta = 300 to 360 degrees

! Each part has (n+1)*(n+2)/2 nodes: 1+2+3+...+(nr_gs+1) = (n_gs+1)*(n_gs+2)/2
  allocate(node_map((nr_gs + 1)*(nr_gs + 2)/2))

 new_sectors : do i = 1, 5

! (1)Generate new nodes and construct a node map (one to one)
!    inode = local node numbering for node_map(:)

  node_map(1) = 1 !Node at the origin

  do k = 2, nr_gs + 1 !Origin to outer boundary
   do j = 1, k        !Right to left

    inode = (k-1)*k/2 + j !Local node number = Right-most node + increment(1,2,..,k)

!   Right side nodes are existing nodes: Left side nodes of the previous sector
    if (j==1) then

     if     (i == 1) then
      node_map(inode) = (k-1)*k/2 + k         !Left-most node of the original sector
     elseif (i == 2) then
      node_map(inode) = nnodes_gs + (k-1)*k/2 !Left-most node of the second sector
     else
      node_map(inode) = nnodes_gs + (nnodes_gs-(nr_gs+1))*(i-2) + (k-1)*k/2
     endif

!   Left side of the last one is the right side of the original sector
    elseif (i==5 .and. j==k) then

      node_map(inode) = (k-1)*k/2 + 1         !Right-most node of the original sector

!   New nodes: Rotate the original nodes by theta = i*pi/3 (i times 60 degrees).
    else

     theta = (pi/three) * real(i,dp)
     nnodes_disk = nnodes_disk + 1
     node_map(inode) = nnodes_disk
     node_disk(nnodes_disk)%x = cos(theta)*node_gs(inode)%x - sin(theta)*node_gs(inode)%y
     node_disk(nnodes_disk)%y = sin(theta)*node_gs(inode)%x + cos(theta)*node_gs(inode)%y

     if (i==1) then

      k1_disk(nnodes_disk) = -k2_gs(inode)
      k2_disk(nnodes_disk) = -k3_gs(inode)
      k3_disk(nnodes_disk) = -k1_gs(inode)

     elseif(i==2) then

      k1_disk(nnodes_disk) =  k3_gs(inode)
      k2_disk(nnodes_disk) =  k1_gs(inode)
      k3_disk(nnodes_disk) =  k2_gs(inode)
 
     elseif(i==3) then
 
      k1_disk(nnodes_disk) = -k1_gs(inode)
      k2_disk(nnodes_disk) = -k2_gs(inode)
      k3_disk(nnodes_disk) = -k3_gs(inode)

     elseif(i==4) then
 
      k1_disk(nnodes_disk) =  k2_gs(inode)
      k2_disk(nnodes_disk) =  k3_gs(inode)
      k3_disk(nnodes_disk) =  k1_gs(inode)

     elseif(i==5) then
 
      k1_disk(nnodes_disk) = -k3_gs(inode)
      k2_disk(nnodes_disk) = -k1_gs(inode)
      k3_disk(nnodes_disk) = -k2_gs(inode)

     endif

!    Keep the record of corner nodes and their neighbors
     if (trim(remove_corner_nodes) == "yes") then
      if (j == k .and. k == nr_gs+1 .and. i < 5) nodes_removed(i+2)       = nnodes_disk
      if (j == k .and. k == nr_gs   .and. i < 5) nodes_removed_nghbr(i+2) = nnodes_disk
     endif

    endif
   end do
  end do

! (2)Generate triangles on the new sector.

  do k = 1, ntrias_gs
   ntrias_disk = ntrias_disk + 1
   tria_disk(ntrias_disk)%v    = node_map(tria_gs(k)%v)
   tria_disk(ntrias_disk)%type = tria_gs(k)%type
  end do

 end do new_sectors

  write(*,*) " nnodes_disk = ", nnodes_disk
  write(*,*) " ntrias_disk = ", ntrias_disk

!--------------------------------------------------------------------------------
! (If requested) Remove the six nodes (or equivalently remove 12 triangles)
! to avoid locally small cells typically caused by smoothing.
!
 if (trim(remove_corner_nodes) == "yes" .and. n_sections == 6) then

 remove_node : do k = 1, 6
    rnode = nodes_removed(k)
  newnode = nodes_removed_nghbr(k)

! (1)Remove the "rnode"

!  Move newnode to the boundary (rnode position)
   node_disk(newnode)%x = node_disk(rnode)%x
   node_disk(newnode)%y = node_disk(rnode)%y
!  Change the last node index to rnode to save/keep the last node
   node_disk(rnode)%x   = node_disk(nnodes_disk)%x
   node_disk(rnode)%y   = node_disk(nnodes_disk)%y
   nodes_changed(k,1)   = nnodes_disk
   nodes_changed(k,2)   = rnode
!  Reduce the dimension: "rnode" information is gone now.
   nnodes_disk = nnodes_disk - 1

! (2)Delete triangles with "rnode" and modify those having "nnodes_disk+1"

   i = 0
  do
   i = i + 1
   if ( i > ntrias_disk) exit

!  Delete if the triangle has "rnode" as its vertex.
   deleted = .false.
   if (tria_disk(i)%v(1) == rnode .or. &
       tria_disk(i)%v(2) == rnode .or. &
       tria_disk(i)%v(3) == rnode        ) then

       tria_disk(i)%v    = tria_disk(ntrias_disk)%v 
       tria_disk(i)%type = tria_disk(ntrias_disk)%type
       ntrias_disk = ntrias_disk - 1
       deleted = .true.
   endif

!  Check if the triangle has the node "nnodes_disk+1" (which no longer exists)
!  as its vertex. If so, overwrite it by rnode.
   if(.not.deleted) then

    if     (tria_disk(i)%v(1) == nnodes_disk+1) then
     tria_disk(i)%v(1) = rnode
    elseif (tria_disk(i)%v(2) == nnodes_disk+1) then
     tria_disk(i)%v(2) = rnode
    elseif (tria_disk(i)%v(3) == nnodes_disk+1) then
     tria_disk(i)%v(3) = rnode
    endif

   endif

!  If the triangle is removed, then shift the index to inspect the new
!  triangle ('i' is again a new triangle coming from ntrias_disk).
   if (deleted) i = i - 1

  end do

 end do remove_node

  write(*,*) " 6 corner nodes have been removed."
  write(*,*) " Updated nnodes_disk = ", nnodes_disk
  write(*,*) " Updated ntrias_disk = ", ntrias_disk

 endif

!--------------------------------------------------------------------------------
! (If requested) Apply smoothing on the disk triangulation

 smooth_it : if (smoothing_type /= 0 .and. n_sections == 6) then

  write(*,*) "Applying smoothing..."

  allocate(dxy(nnodes_disk,2))
  smoothing_factor = 0.05_dp
     max_smoothing = 1000

 smoothing : do i_smoothing = 1, max_smoothing

  dxy = zero !Initialize the changes

! Accumulate the changes by looping over triangles
  do i = 1, ntrias_disk

   v = tria_disk(i)%v

!  x-coordinate
    dxy(v(1),1) = dxy(v(1),1) + ( node_disk(v(2))%x - node_disk(v(1))%x )
    dxy(v(1),1) = dxy(v(1),1) + ( node_disk(v(3))%x - node_disk(v(1))%x )

    dxy(v(2),1) = dxy(v(2),1) + ( node_disk(v(1))%x - node_disk(v(2))%x )
    dxy(v(2),1) = dxy(v(2),1) + ( node_disk(v(3))%x - node_disk(v(2))%x )

    dxy(v(3),1) = dxy(v(3),1) + ( node_disk(v(1))%x - node_disk(v(3))%x )
    dxy(v(3),1) = dxy(v(3),1) + ( node_disk(v(2))%x - node_disk(v(3))%x )


!  y-coordinate
    dxy(v(1),2) = dxy(v(1),2) + ( node_disk(v(2))%y - node_disk(v(1))%y )
    dxy(v(1),2) = dxy(v(1),2) + ( node_disk(v(3))%y - node_disk(v(1))%y )

    dxy(v(2),2) = dxy(v(2),2) + ( node_disk(v(1))%y - node_disk(v(2))%y )
    dxy(v(2),2) = dxy(v(2),2) + ( node_disk(v(3))%y - node_disk(v(2))%y )


    dxy(v(3),2) = dxy(v(3),2) + ( node_disk(v(1))%y - node_disk(v(3))%y )
    dxy(v(3),2) = dxy(v(3),2) + ( node_disk(v(2))%y - node_disk(v(3))%y )

  end do

! Make changes to each node except the boundary nodes.
  dxy_norm = -1000000.0_dp
  do i=1, nnodes_disk

   if ( abs( sqrt(node_disk(i)%x**2 + node_disk(i)%y**2) - r_gs ) < 1.0e-14 ) then
   else

!  Constrained smoothing: skip nodes in the sector boundaries.
   if (smoothing_type == 2) then

    if ( abs(atan2(node_disk(i)%y, node_disk(i)%x) - pi/three)     < 1.0e-13 ) cycle
    if ( abs(atan2(node_disk(i)%y, node_disk(i)%x) + pi/three)     < 1.0e-13 ) cycle
    if ( abs(atan2(node_disk(i)%y, node_disk(i)%x) - two*pi/three) < 1.0e-13 ) cycle
    if ( abs(atan2(node_disk(i)%y, node_disk(i)%x) + two*pi/three) < 1.0e-13 ) cycle
    if ( abs(node_disk(i)%y) < 1.0e-13 ) cycle

   endif

     node_disk(i)%x = node_disk(i)%x + smoothing_factor * dxy(i,1)
     node_disk(i)%y = node_disk(i)%y + smoothing_factor * dxy(i,2)

!    L_inf norm of changes scaled by the typical mesh spacing, dr_gs.
     dxy_norm = max( dxy_norm, abs(dxy(i,1)**2 + dxy(i,2)**2)/dr_gs )

   endif

 end do

! Exit if converged
  if ( dxy_norm < 1.0e-04) then
   write(*,*) " Smoothing converged at ", i_smoothing
   exit smoothing
  elseif (i_smoothing == max_smoothing) then
   write(*,*) " Smoothing didn't converge... ", "  dxy_norm = ", dxy_norm
  endif

 end do smoothing

 deallocate(dxy)

 endif smooth_it

 deallocate(k1_gs,k2_gs,k3_gs,node_gs,tria_gs,node_map)

!--------------------------------------------------------------------------------

! Now, at this point, we have a triangulation of a disk defined by
! tria_disk and node_disk.
! Number of triangles = ntrias_disk
! Number of nodes     = nnodes_disk

!*******************************************************************************
! Write a Tecplot file for the triangulated disk.
!******************************************************************************
 debug_mode_02 : if (debug_mode) then
 open(unit=2, file=filename_disk, status="unknown", iostat=os)

  write(2,*) 'TITLE = "GRID"'
  write(2,*) 'VARIABLES = "x","y","z","k1","k2","k3","kd"'
  write(2,*) 'ZONE  N=', nnodes_disk,',E=', ntrias_disk,' , ET=triangle, F=FEPOINT'

! Nodes
  do i = 1, nnodes_disk

   sk12p = ( 1 + int_sign( k1_disk(i)*k2_disk(i) ) )/2
   sk23p = ( 1 + int_sign( k2_disk(i)*k3_disk(i) ) )/2
   sk31p = ( 1 + int_sign( k3_disk(i)*k1_disk(i) ) )/2

   sk12m = ( 1 - int_sign( k1_disk(i)*k2_disk(i) ) )/2
   sk23m = ( 1 - int_sign( k2_disk(i)*k3_disk(i) ) )/2
   sk31m = ( 1 - int_sign( k3_disk(i)*k1_disk(i) ) )/2

   if     (sk12p > 0) then

     sk23p = 0
     sk23m = 1

     sk31p = 0
     sk31m = 1

   elseif (sk23p > 0) then

     sk12p = 0
     sk12m = 1

     sk31p = 0
     sk31m = 1

   elseif (sk31p > 0) then

     sk12p = 0
     sk12m = 1

     sk23p = 0
     sk23m = 1

   endif

    write(2,'(3ES20.10,4i10)') node_disk(i)%x, node_disk(i)%y, 0.0, k1_disk(i),k2_disk(i),k3_disk(i), &
    sk31m*sk23m*sk12p*abs(k3_disk(i)) + &
    sk31m*sk12m*sk23p*abs(k1_disk(i)) + &
    sk23m*sk12m*sk31p*abs(k2_disk(i))

  end do

! Triangles
  do i = 1, ntrias_disk
   write(2,'(3I10)') tria_disk(i)%v(1), tria_disk(i)%v(2), tria_disk(i)%v(3)
  end do

 close(2)

 write(*,*) "Tecplot file has been written: ", filename_disk
 endif debug_mode_02
!*******************************************************************************

!*******************************************************************************
! 3. Map the disk triangulation onto the hemisphere surface.
!
! OK, now, let's gently attach the disk onto the hemisphere.
! This is done locally at each node, just by using the node_disk data.
! It doesn't matter how they are ordered.
! Connectivity data are unchanged.
!*******************************************************************************

  write(*,*) "***********************************************************"
  write(*,*) " 3. Place the disk onto the hemisphere "
  write(*,*) "***********************************************************"

  nnodes_circum = n_sections*nr_gs

  if (n_sections == 1) nnodes_circum = nnodes_circum + 1

  nnodes = nnodes_disk
  allocate(node_body_hs(nnodes))
  allocate(k1_body_hs(nnodes))
  allocate(k2_body_hs(nnodes))
  allocate(k3_body_hs(nnodes))
  allocate(k4_body_hs(nnodes))

   nnodes_hs = 0

  do i = 1, nnodes_disk

   k1_body_hs(i) = k1_disk(i) ! Flip the sign JFC
   k2_body_hs(i) = k2_disk(i) ! Flip the sign JFC
   k3_body_hs(i) = k3_disk(i) ! Flip the sign JFC
   k4_body_hs(i) =  0          ! <- On the hemisphere

! Push the node onto the circle located at z=0.0.
   s = sqrt(node_disk(i)%x**2 + node_disk(i)%y**2)
   if (i==1) then
    xp = zero
    yp = zero
   else
    xp = node_disk(i)%x/s * Rd !Extend it to the circle of radius Rd
    yp = node_disk(i)%y/s * Rd !Extend it to the circle of radius Rd
   endif

!    zp = 0.0  The circle is located at x = 0.0.

! Now, the node (xp,yp,zp) is located on the perimeter of the circle at z=0.0.
! Rotate the node onto the sphere.
            theta = s/Rd
           nnodes_hs = nnodes_hs + 1
   node_body_hs(nnodes_hs)%x = xp*sin(theta)
   node_body_hs(nnodes_hs)%y = yp*sin(theta)
   node_body_hs(nnodes_hs)%z = Rd*cos(theta)


! Rotate the hemisphere by 30 degrees.
!  do i = 1, nnodes_hs
   theta = pi/six
   xc = node_body_hs(i)%x
   yc = node_body_hs(i)%y
   node_body_hs(i)%x = cos(theta)*xc - sin(theta)*yc
   node_body_hs(i)%y = sin(theta)*xc + cos(theta)*yc
!  end do


!  Surface normal direction along which we go up to generate prismatic elements.
   node_body_hs(nnodes_hs)%nx = node_body_hs(nnodes_hs)%x
   node_body_hs(nnodes_hs)%ny = node_body_hs(nnodes_hs)%y
   node_body_hs(nnodes_hs)%nz = node_body_hs(nnodes_hs)%z

!  Make it the unit vector (well, probably already a unit vector, though...)
   sf = sqrt(node_body_hs(nnodes_hs)%nx**2 + node_body_hs(nnodes_hs)%ny**2 + node_body_hs(nnodes_hs)%nz**2)
   node_body_hs(nnodes_hs)%nx = node_body_hs(nnodes_hs)%nx / sf
   node_body_hs(nnodes_hs)%ny = node_body_hs(nnodes_hs)%ny / sf
   node_body_hs(nnodes_hs)%nz = node_body_hs(nnodes_hs)%nz / sf

  end do

 deallocate(k1_disk,k2_disk,k3_disk,node_disk)
 write(*,*) " nnodes_hs, nnodes_disk = ", nnodes_hs, nnodes_disk

!*******************************************************************************
! Write a Tecplot file for the triangulated hemisphere surface.
!******************************************************************************
 debug_mode_03 : if (debug_mode) then
 open(unit=3, file=filename_surface, status="unknown", iostat=os)

  write(3,*) 'TITLE = "GRID"'
  write(3,*) 'VARIABLES = "x","y","z","k1+k2"'
  write(3,*) 'ZONE  N=', nnodes_hs,',E=', ntrias_disk,' , ET=triangle, F=FEPOINT'

! Nodes
  do i = 1, nnodes_hs

    write(3,'(3ES20.10,i10)') node_body_hs(i)%x,  node_body_hs(i)%y, node_body_hs(i)%z, &
                          k1_body_hs(i)+ k2_body_hs(i)
  end do

! Triangles
  do i = 1, ntrias_disk
   write(3,'(3I10)') tria_disk(i)%v(1), tria_disk(i)%v(2), tria_disk(i)%v(3)
  end do

 close(3)

 write(*,*) "Tecplot file has been written: ", filename_surface
 endif debug_mode_03
!*******************************************************************************

!*******************************************************************************
! 4. Create some temporary data, and modify some data.
!
!    At this point, the surface grid is still a hemisphere.
!
!                  z
!                  ^
!        o  o      |
!     o        o   |
!    o__________o  -------> x
!               
!
!*******************************************************************************
!-------------------------------------------------------------------------------
! Construct a 2D node list around the circle of the bottom of the hemisphere.
!
!  nodes_circum(k), k=1,nnodes_circum
!
!         <-
!        o  o    
!     o        o     y
!    o          o    ^
!    o          o    |
!     o        o     |
!        o  o        -------> x
!         ->

  allocate(nodes_circum(nnodes_circum+1))
  write(*,*) " nnodes around the circle of the hemisphere = ", nnodes_circum

  nnodes_circum = 0

! Circum of generating sector: Counter-clockwise ordering.
!
  do k = 1, nr_gs+1
   nnodes_circum = nnodes_circum + 1
   nodes_circum(nnodes_circum) = nr_gs*(nr_gs+1)/2 + k
  end do

! Circum of other sectors
  do i = 1, 4
   do k = 1, nr_gs
    nnodes_circum = nnodes_circum + 1
    nodes_circum(nnodes_circum) = nnodes_gs + (i-1)*(nr_gs+1)*nr_gs/2 + (nr_gs-1)*nr_gs/2 + k
   end do
  end do

! i == 5
   do k = 1, nr_gs-1
    nnodes_circum = nnodes_circum + 1
    nodes_circum(nnodes_circum) = nnodes_gs + (5-1)*(nr_gs+1)*nr_gs/2 + (nr_gs-2)*(nr_gs-1)/2 + k
   end do

  nodes_circum(nnodes_circum + 1) = nodes_circum(1)

! Fix the node number if 6 nodes have been removed.

  if (trim(remove_corner_nodes) == "yes" .and. n_sections == 6) then

   do i = 1, nnodes_circum + 1

    do k = 1, 6

     if (nodes_circum(i) == nodes_removed(k)) then

      nodes_circum(i) = nodes_removed_nghbr(k)

     elseif (nodes_circum(i) == nodes_changed(k,1)) then

      nodes_circum(i) = nodes_changed(k,2)

     endif

    end do

   end do

  endif

!-------------------------------------------------------------------------------
! Generate new array for triangles.

  ntrias_hs = ntrias_disk
  allocate(tria_hs(ntrias_hs))

! Copy the data for triangles.
! Reverse the node-ordering.
! Triangles are now ordered in counter-clockwise manner,
! so that it points towards inside the domain.
!
!                            1
!                            o
!  Interior domain         . .
!                     <---.- .
!                        .   .
!                       o----o
!                      2     3

  do i = 1, ntrias_hs

   node1 = tria_disk(i)%v(1)
   node2 = tria_disk(i)%v(2)
   node3 = tria_disk(i)%v(3)

    tria_hs(i)%v(1) = node3
    tria_hs(i)%v(2) = node2
    tria_hs(i)%v(3) = node1

   tria_hs(i)%type = tria_disk(i)%type

  end do

  deallocate(tria_disk)

!*******************************************************************************
! Generate a spherical surface grid by duplicating and combining the hemisphere
! surface grids. 
!
!                  z
!                  ^
!        o  o      |
!     o        o   |
!    o__________o  -------> x
!    o          o
!     o        o
!        o  o
!               
!
!*******************************************************************************

  nnodes = 2*nnodes_hs - nnodes_circum
  allocate(node_body(nnodes))
  allocate(k1_body(nnodes))
  allocate(k2_body(nnodes))
  allocate(k3_body(nnodes))
  allocate(k4_body(nnodes))

! Copy the hemisphere surface data.

  do i = 1, nnodes_hs

   node_body(i)%x   = node_body_hs(i)%x
   node_body(i)%y   = node_body_hs(i)%y
   node_body(i)%z   = node_body_hs(i)%z

   node_body(i)%nx  = node_body_hs(i)%nx
   node_body(i)%ny  = node_body_hs(i)%ny
   node_body(i)%nz  = node_body_hs(i)%nz

   node_body(i)%nx2 = node_body_hs(i)%nx2
   node_body(i)%ny2 = node_body_hs(i)%ny2
   node_body(i)%nz2 = node_body_hs(i)%nz2

   k1_body(i) = k1_body_hs(i)
   k2_body(i) = k2_body_hs(i)
   k3_body(i) = k3_body_hs(i)
   k4_body(i) = k4_body_hs(i)

  end do

  ntrias = 2*ntrias_hs
  allocate(tria(ntrias))

  do i = 1, ntrias_hs
   tria(i)%v(1) = tria_hs(i)%v(1)
   tria(i)%v(2) = tria_hs(i)%v(2)
   tria(i)%v(3) = tria_hs(i)%v(3)
   tria(i)%type = tria_hs(i)%type
  end do

! Generate the other half:
!   Flip the hemisphere wrt z-axis and add it to
!   the original hemisphere.

  do i = 1, nnodes_hs
   node_body_hs( i )%circum_node = .false.
  end do

  do i = 1, nnodes_circum
   node_body_hs( nodes_circum(i) )%circum_node = .true.
  end do

  allocate(node_map(nnodes_hs))

  nnodes = nnodes_hs
 
  do i = 1, nnodes_hs

   if ( node_body_hs(i)%circum_node ) then
    node_map(i)           = i
    cycle
   endif

                  nnodes =   nnodes + 1
   node_map(i)           =   nnodes

   node_body(nnodes)%x   =   node_body_hs(i)%x
   node_body(nnodes)%y   =   node_body_hs(i)%y
   node_body(nnodes)%z   = - node_body_hs(i)%z ! Flip

   node_body(nnodes)%nx  =   node_body_hs(i)%nx
   node_body(nnodes)%ny  =   node_body_hs(i)%ny
   node_body(nnodes)%nz  = - node_body_hs(i)%nz ! Flip

   node_body(nnodes)%nx2 =   node_body_hs(i)%nx2
   node_body(nnodes)%ny2 =   node_body_hs(i)%ny2
   node_body(nnodes)%nz2 = - node_body_hs(i)%nz2 ! Flip

   k1_body(nnodes) = k1_body_hs(i)
   k2_body(nnodes) = k2_body_hs(i)
   k3_body(nnodes) = k3_body_hs(i)
   k4_body(nnodes) = k4_body_hs(i)

  end do

  ntrias = ntrias_hs

  do i = 1, ntrias_hs

         ntrias = ntrias + 1

   tria(ntrias)%v(1) = node_map( tria_hs(i)%v(1) )
   tria(ntrias)%v(2) = node_map( tria_hs(i)%v(3) )
   tria(ntrias)%v(3) = node_map( tria_hs(i)%v(2) )
   tria(ntrias)%type = tria_hs(i)%type

  end do

! Rotate the sphere wrt y-axis: 90 degrees anticlockwise in (x,z) plane around the origin.

  do i = 1, nnodes

   theta = pi/two
   xc = node_body(i)%x
   zc = node_body(i)%z
   node_body(i)%x = cos(theta)*xc - sin(theta)*zc
   node_body(i)%z = sin(theta)*xc + cos(theta)*zc

!  Surface normal direction along which we go up to generate prismatic elements.
   node_body(i)%nx = node_body(i)%x
   node_body(i)%ny = node_body(i)%y
   node_body(i)%nz = node_body(i)%z

!  Make it the unit vector (well, probably already a unit vector, though...)
   sf = sqrt(node_body(i)%nx**2 + node_body(i)%ny**2 + node_body(i)%nz**2)
   node_body(i)%nx = node_body(i)%nx / sf
   node_body(i)%ny = node_body(i)%ny / sf
   node_body(i)%nz = node_body(i)%nz / sf

  end do

  deallocate(node_body_hs,tria_hs,k1_body_hs,k2_body_hs,k3_body_hs,k4_body_hs)

!--------------------------------------------------------------------------------
! (If requested) Apply smoothing on the surface triangulation

 smooth_it_surface : if (smoothing_type /= 0 .and. n_sections == 6) then

  write(*,*)
  write(*,*) "  Applying smoothing on the sphere surface..."

  allocate(dxyz(nnodes,3))
  smoothing_factor = 0.05_dp
     max_smoothing = 1000

 smoothing_surface : do i_smoothing = 1, max_smoothing

  dxyz = zero !Initialize the changes

! Accumulate the changes by looping over triangles
  do i = 1, ntrias

   v = tria(i)%v

!  x-coordinate
    dxyz(v(1),1) = dxyz(v(1),1) + ( node_body(v(2))%x - node_body(v(1))%x )
    dxyz(v(1),1) = dxyz(v(1),1) + ( node_body(v(3))%x - node_body(v(1))%x )

    dxyz(v(2),1) = dxyz(v(2),1) + ( node_body(v(1))%x - node_body(v(2))%x )
    dxyz(v(2),1) = dxyz(v(2),1) + ( node_body(v(3))%x - node_body(v(2))%x )

    dxyz(v(3),1) = dxyz(v(3),1) + ( node_body(v(1))%x - node_body(v(3))%x )
    dxyz(v(3),1) = dxyz(v(3),1) + ( node_body(v(2))%x - node_body(v(3))%x )


!  y-coordinate
    dxyz(v(1),2) = dxyz(v(1),2) + ( node_body(v(2))%y - node_body(v(1))%y )
    dxyz(v(1),2) = dxyz(v(1),2) + ( node_body(v(3))%y - node_body(v(1))%y )

    dxyz(v(2),2) = dxyz(v(2),2) + ( node_body(v(1))%y - node_body(v(2))%y )
    dxyz(v(2),2) = dxyz(v(2),2) + ( node_body(v(3))%y - node_body(v(2))%y )

    dxyz(v(3),2) = dxyz(v(3),2) + ( node_body(v(1))%y - node_body(v(3))%y )
    dxyz(v(3),2) = dxyz(v(3),2) + ( node_body(v(2))%y - node_body(v(3))%y )


!  z-coordinate
    dxyz(v(1),3) = dxyz(v(1),3) + ( node_body(v(2))%z - node_body(v(1))%z )
    dxyz(v(1),3) = dxyz(v(1),3) + ( node_body(v(3))%z - node_body(v(1))%z )

    dxyz(v(2),3) = dxyz(v(2),3) + ( node_body(v(1))%z - node_body(v(2))%z )
    dxyz(v(2),3) = dxyz(v(2),3) + ( node_body(v(3))%z - node_body(v(2))%z )

    dxyz(v(3),3) = dxyz(v(3),3) + ( node_body(v(1))%z - node_body(v(3))%z )
    dxyz(v(3),3) = dxyz(v(3),3) + ( node_body(v(2))%z - node_body(v(3))%z )

  end do

! Make changes to each node except the boundary nodes.
  dxyz_norm = -1000000.0_dp
  do i = 1, nnodes

     xc = node_body(i)%x
     yc = node_body(i)%y
     zc = node_body(i)%z

     node_body(i)%x = node_body(i)%x + smoothing_factor * dxyz(i,1)
     node_body(i)%y = node_body(i)%y + smoothing_factor * dxyz(i,2)
     node_body(i)%z = node_body(i)%z + smoothing_factor * dxyz(i,3)

   !Project back onto the surface.
     sf = sqrt(node_body(i)%x**2 + node_body(i)%y**2 + node_body(i)%z**2)
     node_body(i)%x = node_body(i)%x/sf * Rd
     node_body(i)%y = node_body(i)%y/sf * Rd
     node_body(i)%z = node_body(i)%z/sf * Rd

!    L_inf norm of changes.
     dxyz_norm = max( dxyz_norm, abs((xc-node_body(i)%x)**2         &
               + (yc-node_body(i)%y)**2 + (zc-node_body(i)%z)**2) )

 end do

! Exit if converged
  if ( dxyz_norm < 1.0e-04) then
   write(*,*) " Smoothing converged at ", i_smoothing
   exit smoothing_surface
  elseif (i_smoothing == max_smoothing) then
   write(*,*) " Smoothing didn't converge... ", "  dxyz_norm = ", dxyz_norm
  endif

 end do smoothing_surface

 deallocate(dxyz)

  write(*,*)

! Re-compute the normals.

  do i = 1, nnodes
   node_body(i)%nx = node_body(i)%x
   node_body(i)%ny = node_body(i)%y
   node_body(i)%nz = node_body(i)%z
   sf = sqrt(node_body(i)%nx**2 + node_body(i)%ny**2 + node_body(i)%nz**2)
   node_body(i)%nx = node_body(i)%nx / sf
   node_body(i)%ny = node_body(i)%ny / sf
   node_body(i)%nz = node_body(i)%nz / sf
  end do

 endif smooth_it_surface


!*******************************************************************************
! Write a Tecplot file for the triangulated sphere surface.
!******************************************************************************
 debug_mode_033 : if (debug_mode) then
 open(unit=3, file=filename_surface_s, status="unknown", iostat=os)

  write(3,*) 'TITLE = "GRID"'
  write(3,*) 'VARIABLES = "x","y","z","k1+k2"'
  write(3,*) 'ZONE  N=', nnodes,',E=', ntrias,' , ET=triangle, F=FEPOINT'

! Nodes
  do i = 1, nnodes

    write(3,'(3ES20.10,i10)') node_body(i)%x,  node_body(i)%y, node_body(i)%z, &
                          k1_body(i)+ k2_body(i)
  end do

! Triangles
  do i = 1, ntrias
   write(3,'(3I10)') tria(i)%v(1), tria(i)%v(2), tria(i)%v(3)
  end do

 close(3)

 write(*,*) "Tecplot file has been written: ", filename_surface
 endif debug_mode_033
!*******************************************************************************

!*******************************************************************************
! 5. We now go up to generate the interior nodes.
!
!                  z
!                  ^
!        o  o      |
!     o        o   -------> x
!    o__________o------------------------------->o   Outer boundary point.
!    o          o  Go up and generate nodes ->
!     o        o
!        o  o
!
! Determine the vertical spacing, vspacing: r(i+1) = vspacing(i) + r(i),
! and generate interior nodes by going up in the direction of surface normal.
!
!   Geometric sequence       Exponential stretching
! o-o--o---o----o-----o-----x-------x----------x--------------x        -----> r
! 1 2  3              nm                                 nr
!   Prismatic layer       Outer tetrahedral region
!
! NOTE: Here, r indicates the direction of surface normal.
! NOTE: The number of nodes in the prismatic layer (nm) is automatically determined
!       for a specified the aspect ratio of the first element on the body = 50
!
! NOTE: The exponential stretching is determined by iterations to equate
!       the spacing of the elements below and above the nm-th node.
!
!*******************************************************************************

  write(*,*) "***********************************************************"
  write(*,*) " 5. Go up to the outer boundary and generate interior nodes"
  write(*,*) "***********************************************************"

!  Estimate of the arc length of an edge along the hemisphere.

    ds = (half*pi*Rd) / real(nr_gs,dp)
    write(*,*) "Mesh spacing along the hemisphere (body)= ", ds

!  Estimate of the arc length of one cell along the outer boundary.

    ds_outer = distance*(ds/Rd)
    write(*,*) "Mesh spacing along the hemisphere (farfirled) = ", ds_outer

! rmax = thickness of the prismatic region (boundary layer):
! If a negative value is given, we set this to be 0.1:
  if (rmax < zero) rmax = 1.0e-01_dp

  write(*,*) " Boundary layer grid thickness = ", rmax

! Compute the parameters for geometric sequence.
! Determine the number of nodes in the r-direction for a given max AR.

  drN = ds !* 0.5_dp !Spacing of the outer most cell in the layer.

  if (drN > rmax) then
   write(*,*) ">>>>> drN > rmax !! Impossible... drN = ds = ", drN
   drN = rmax *0.5_dp
   write(*,*) ">>>>> Modified drN = rmax *0.5_dp = ", drN
  endif

!---------------------------------------------------------------------
! dr1 is the first vertical spacing off the wall.

  set_dr1 : if ( target_reynolds_number > zero) then

   write(*,*) ">>>>>>>>>>>>>>> dr1 from a target Re and y_plus:"
    cf = 0.026_dp/target_reynolds_number**(1.0_dp/7.0_dp)
   dr1 = ( sqrt(2.0_dp/cf)/target_reynolds_number) * target_yplus
   write(*,*) "--- dr1 determined by the target Re and y_plus:"
   write(*,*) " target_reynolds_number = ", target_reynolds_number
   write(*,*) " target_yplus = ", target_yplus
   write(*,*) " --> dr1 = ", dr1
   write(*,*) "      AR = ", ds/dr1
   write(*,*) ">>>>>>>>>>>>>>> "

  else

   dr1 = target_yplus

  endif set_dr1
!---------------------------------------------------------------------

   if (dr1 > rmax) then
    dr1 = 0.1_dp * ds
    write(*,*) " dr1 too large: >> rmax = ", rmax
    write(*,*) " --> Adjusted: dr1 = ", dr1 
   endif

   gr = (rmax-dr1)/(rmax-drN)
   write(*,*) "      gr = ", gr
   nm = ceiling( log(gr*drN/dr1)/log(gr) )
   gr = (drN/dr1)**( one/(real(nm,dp)-one) )

   write(*,*)
   write(*,*) " NOTE: dr1 is a non-dimensionalized length. "
   write(*,*) "       Reference length is the diameter of the hemisphere (=1.0 in the grid)."
   write(*,*)

   write(*,*) "ds (mesh spacing along the hemisphere) = ", ds
   write(*,*) "rmax (the end of the boundary layer)      = ", rmax
   write(*,*) "dr1 = ", dr1
   write(*,*) "drN = ", drN
   write(*,*) "gr = ", gr
   write(*,*) "dr1*gr**(nm-1) = ", dr1*gr**(nm-1)

!  At this point, nm is the number of elements within the layer.
!  Add 1 to make it the number of nodes.
   nm = nm + 1 ! where change the type of cell from prism to tetra for a mixed grid.
   write(*,*) "nm (nodes) = ", nm

!  Allocate the vertical spacing array
   allocate(vspacing_nm(nm-1))

  do i = 1, nm-1

   if (i == 1) then

    vspacing_nm(i) = dr1

!  Geometric sequence for the cell spacing inside the layer.
   elseif (i < nm) then

    vspacing_nm(i) = dr1*gr**(i-1)

   endif

  end do

  write(*,*)

! NOTE: The last node in the prismatic layer is nm-th node.
!       The spacing array vspacing() goes from 1 to nm-1.

!  Take a guess on the number of nodes in the outer tet region.

     nr = nm + int(ceiling( (distance - rmax) / ds )/25.0_dp)

   write(*,*) "Initial value of nr = ", nr

   rs_nmm1 = zero
  do i = 1, nm-2
   rs_nmm1 = rs_nmm1 + vspacing_nm(i)
  end do

! Uniform spacing coordinate in the outer region for now.
   rs_nm   = rs_nmm1 + vspacing_nm(nm-1)
  write(*,*) "  rs_nmm1 = ", rs_nmm1
  write(*,*) "  rs_nm   = ", rs_nm

  write(*,*)
  write(*,*) ">>>>>>>>>>>>>>> Start iteration for nr."
!-------------------------------------------------------------------
 determine_nr_itaratively : do i_nr = 1, 500
!-------------------------------------------------------------------

  write(*,*)

   dr_outer = (distance - rs_nm)/real(nr-nm,dp)
   rs_nmp1 = rs_nm   + dr_outer
  write(*,*) "  rs_nmp1 = ", rs_nmp1

!-------------------------------------------------------------------
! Compute the stretching factor for a smooth transition from
! the boundary layer to outer region.
!  Target spacing: 1.3 times the spacing below for a smooth transition
    s = (rs_nm-rs_nmm1) * 3.0_dp
!  Initial value
   sf = 2.0_dp
!  Determine sf iteratively
   do k = 1, 500
    drnmp1 = (one-exp(sf*(rs_nmp1-rs_nm)/(distance-rs_nm)))/(one-exp(sf)) * (distance-Rd)
    if (abs(drnmp1 - s)/s < 1.0e-02_dp) exit
    if (drnmp1 > s ) then
      sf = sf + 0.01_dp
    else
      sf = sf - 0.01_dp
    endif
   end do
!-------------------------------------------------------------------

  write(*,*) "Determined stretching factor = ",sf, " at iterations = ", k

  rs_nr   = distance-Rd
  rs_nrm1 = rs_nr - dr_outer
  rs_nrm1 = rs_nm + (one-exp(sf*(rs_nrm1-rs_nm)/(distance-rs_nm)))/(one-exp(sf)) * (distance-Rd)

  write(*,*) "  rs_nr   = ", rs_nr
  write(*,*) "  rs_nrm1 = ", rs_nrm1

  dr_spacing = rs_nr - rs_nrm1
  write(*,*) "  dr_spacing    = ", dr_spacing 
  err_dr_spacing = (dr_spacing-ds_outer)/ds_outer
  write(*,*) "  Difference(dr_spacing,ds_outer)=", err_dr_spacing

  if (abs(err_dr_spacing) < 0.30_dp) then
   write(*,*)
   write(*,*) " Iteration converged for nr: final nr = ", nr
   exit determine_nr_itaratively
  endif

 !If not satisfactory, adjust nr and try again.
  if (err_dr_spacing < 0.0_dp) then
   nr = nr - 2
  else
   nr = nr + 2
  endif

!-------------------------------------------------------------------
  end do determine_nr_itaratively
!-------------------------------------------------------------------

  write(*,*)
  write(*,*) ">>>>>>>>>>>>>>> End of iteration for nr."

! Create a vertical coordinate array.
! Note: rs=0        on the hemisphere surface.
!       rs=distance at the outer boudnary.
  allocate(rs(nr))
  rs(1) = zero

! Compute the vertical coordinate based on the spacing generated by
! the geometric sequence as above.
  do i = 2, nm
   rs(i) = vspacing_nm(i-1) + rs(i-1)
  end do

  allocate(vspacing(nr))
  do i = 1, nm-1
   vspacing(i) = vspacing_nm(i)
  end do
   vspacing(nm) = vspacing_nm(nm-1)

  deallocate(vspacing_nm)

! Uniform spacing coordinate:
   dr_outer = (distance - rs(nm))/real(nr-nm,dp)
  do i = nm+1, nr
   rs(i) = rs(nm) + real(i-nm,dp)*dr_outer
  end do

  do i = nm+1, nr
   rs(i) = rs(nm) + (one-exp(sf*(rs(i)-rs_nm)/(distance-rs_nm)))/(one-exp(sf)) * (distance-Rd)
  end do
  !Adjust the last one, so that we have the outer boundary precisely as specified.
   rs(nr) = distance-Rd

! Re-compute the spacing modified by the stretching.
  do i = nm+1, nr
   vspacing(i-1) = rs(i) - rs(i-1)
  end do

! Generate nodes by going up along the surface normal from each node on the body.

  nnodes_body = nnodes
       nnodes = nnodes*nr
  allocate(node(nnodes))
  allocate(node_above(nr*nnodes_body))
  allocate(k1(nnodes))
  allocate(k2(nnodes))
  allocate(k3(nnodes))
  allocate(k4(nnodes))

  write(*,*) " nnodes_body =", nnodes_body 
  do i = 1, nnodes_body
   node(i)%x  = node_body(i)%x
   node(i)%y  = node_body(i)%y
   node(i)%z  = node_body(i)%z
   node(i)%nx = node_body(i)%nx
   node(i)%ny = node_body(i)%ny
   node(i)%nz = node_body(i)%nz
   k1(i)  = k1_body(i)
   k2(i)  = k2_body(i)
   k3(i)  = k3_body(i)
   k4(i)  = k4_body(i)
  end do

  deallocate(node_body)

! We now go up and generate interior nodes

  nnodes = nnodes_body
  node_on_body : do i = 1, nnodes_body

!  Second node
   nnodes = nnodes + 1
   node(nnodes)%x = vspacing(1)*node(i)%nx + node(i)%x
   node(nnodes)%y = vspacing(1)*node(i)%ny + node(i)%y
   node(nnodes)%z = vspacing(1)*node(i)%nz + node(i)%z
   xp = node(nnodes)%x
   yp = node(nnodes)%y
   zp = node(nnodes)%z
   node_above(i) = nnodes

   k1(nnodes)  = k1_body(i)
   k2(nnodes)  = k2_body(i)
   k3(nnodes)  = k3_body(i)
   k4(nnodes)  = 1

!  Third node to nm-th node
!  Nodes in the prismatic layer along the surface normal. 
   do k = 3, nm

    nnodes = nnodes + 1
    node(nnodes)%x = vspacing(k-1)*node(i)%nx + xp
    node(nnodes)%y = vspacing(k-1)*node(i)%ny + yp
    node(nnodes)%z = vspacing(k-1)*node(i)%nz + zp
    node_above(nnodes-1) = nnodes

    xp = node(nnodes)%x
    yp = node(nnodes)%y
    zp = node(nnodes)%z

    k1(nnodes)  = k1_body(i)
    k2(nnodes)  = k2_body(i)
    k3(nnodes)  = k3_body(i)
    k4(nnodes)  = k-1

   end do

!  Outer region
   do k = nm+1, nr

    dirx = node(i)%nx
    diry = node(i)%ny
    dirz = node(i)%nz

!   Go up in the normal direction and generate a new interior node.
    nnodes = nnodes + 1
    node(nnodes)%x = vspacing(k-1)*dirx + xp
    node(nnodes)%y = vspacing(k-1)*diry + yp
    node(nnodes)%z = vspacing(k-1)*dirz + zp
    node_above(nnodes-1) = nnodes

    xp = node(nnodes)%x
    yp = node(nnodes)%y
    zp = node(nnodes)%z

    k1(nnodes)  = k1_body(i)
    k2(nnodes)  = k2_body(i)
    k3(nnodes)  = k3_body(i)
    k4(nnodes)  = k-1

   end do

  end do node_on_body

 deallocate(k1_body,k2_body,k3_body,k4_body)

! At this point, all nodes have been generated.
! It is now a matter of how to connect them (i.e., type of elements).

!*******************************************************************************
! 6. Write a map file: boundary marks
!
! There are four boundary parts:
!
! 1. Hemishpere
! 2. Symmetry wrt z-axis
! 3. Inflow boundary
! 4. Outer boundary
!
! NOTE: The boundary condition numbers (e.g., 4000) are specific to a solver.
!       Appropriate number needs to be assigned for your solver.
!
!*******************************************************************************

  write(*,*) "***********************************************************"
  write(*,*) " 6. Write a boudnary info file"
  write(*,*) "***********************************************************"

  open(unit=16, file=filename_mapbc, status="unknown", iostat=os)

   write(16,'(a57)') "3       !Number of boundary parts (boundary conditions)"
   write(16,'(a21)') "1, 4000 !Viscous wall"
   write(16,'(a15)') "2, 5050 !Inflow"
   write(16,'(a16)') "3, 5051 !Outflow"

  close(16)

  write(*,*) " Boundary info file written:", filename_mapbc

!*******************************************************************************
! 7. Generate elements
!
!    Nodes have already defined.
!    Different grids will be generated by connecting the existing nodes.
!
!    igrid_type = Choice of element type:
!
!      1 = Prismatic grid
!      2 = Tetrahedral grid
!      3 = Mixed grid (Prism/Tet)
!
!*******************************************************************************
!*******************************************************************************
!* (1). Generate a prismatic grid
!*******************************************************************************
 if (igrid_type == 1) then

  write(*,*) "***********************************************************"
  write(*,*) " 7. Generate a prismatic grid"
  write(*,*) "***********************************************************"

 call prismatic_grid

 write(*,*) " Prismatic elements generated."

!*******************************************************************************
!* (2). Generate a tet grid
!*******************************************************************************
 elseif ( igrid_type == 2 ) then

  write(*,*) "***********************************************************"
  write(*,*) " 7. Generate a tetrahedral grid"
  write(*,*) "***********************************************************"

  call tet_grid

 write(*,*) " Tetrahedral elements generated."

!*******************************************************************************
!* (3). Generate a mixed grid  (Prism/Tet)
!*******************************************************************************
 elseif ( igrid_type == 3 ) then

  write(*,*) "***********************************************************"
  write(*,*) " 7. Generate a mixed grid"
  write(*,*) "***********************************************************"

 call mixed_grid

 write(*,*) " Mixed elements generated."

!*******************************************************************************
!* (4). Error
!*******************************************************************************
 else

  write(*,*) "Invalid input: igrid_type = ", igrid_type
  write(*,*) "               igrid_type must be 1, 2, or 3. Try again."
  stop

 endif

!*******************************************************************************
! 9. Generate two line files:
!
!    (1) hemisphere.lines_fmt     - Lines are within the viscous layer.
!    (2) hemisphere.lines_fmt_all - Lines go up to the outer boundary.
!
!*******************************************************************************

  write(*,*) "***********************************************************"
  write(*,*) " 9. Generate lin-info files"
  write(*,*) "***********************************************************"

  if( generate_line_file ) then

! Lines within the viscous layer
  if (line_extent==1) then

  write(*,*)
  write(*,*) "Writing a regular line-info file...."
  call line_info_file(nm)

! Lines go all the way up to the outer boundary (all nodes are in lines)
  elseif (line_extent==2) then

  write(*,*)
  write(*,*) "Writing an all-node line-info file...."
  call line_info_file_all(nr)

  else

  write(*,*) " Invalid input for line_extent:", line_extent
  write(*,*) " line_extent must be 1 or 2..."
  write(*,*) " Try again. Stop."
  stop

  endif
  else
   write(*,*) "Skipping write of line-info files"
   write(*,*) "To enable it, set generate_line_file=.true. with line_extent = 1 or 2."
   write(*,*)
  endif

!*******************************************************************************
! 10. Write a Tecplot file for boundaries.
!*******************************************************************************

  write(*,*) "***********************************************************"
  write(*,*) " 10. Generate a tecplot file for viewing boundaries"
  write(*,*) "***********************************************************"

  if (debug_mode .or. generate_tec_file_b) then

   write(*,*) "Writing a Tecplot file for boundaries...."
   call write_tecplot_boundary_file
   write(*,*) "Tecplot file for boundaries written : ", filename_tecplot_b

  else

   write(*,*) "Skipping write of Teplot boundary file."
   write(*,*) "To enable it, set generate_tec_file_b = .true."
   write(*,*)

  endif

!*******************************************************************************
! 11. Write a k-file.
!*******************************************************************************

  write(*,*) "***********************************************************"
  write(*,*) " 11. Generate a k file"
  write(*,*) "***********************************************************"

  if (generate_k_file) then
    write(*,*) "Writing a k file...."
    call write_k_file
    write(*,*) "K file written : ", filename_k
  else
    write(*,*) "Skipping write of the k file : "
    write(*,*) "To enable it, set generate_k_file = .true."
    write(*,*)
  endif

!*******************************************************************************
! 12. Write a ugrid file.
!*******************************************************************************

  write(*,*) "***********************************************************"
  write(*,*) " 12. Generate a .ugrid file for running a flow solver"
  write(*,*) "***********************************************************"

  write(*,*) "Writing a ugrid file...."
  call write_ugrid_file
  write(*,*) "UGRID file written : ", filename_ugrid

!*******************************************************************************
! 13. Write a Tecplot file for the volume grid.
!*******************************************************************************

   write(*,*) "***********************************************************"
   write(*,*) " 13. Generate a tecplot file for viewing a volume grid"
   write(*,*) "***********************************************************"

  if (generate_tec_file_v) then

   write(*,*) "Writing a Tecplot file for volumes...."
   call write_tecplot_volume_file
   write(*,*) "Tecplot file for volume grid written : ", filename_tecplot_v

  else
   write(*,*) "Skipping write of Teplot volume file."
   write(*,*) "To enable it, set generate_tec_file_v = .true."
  endif

 write(*,*)
 write(*,*) "Congratulations!"
 write(*,*) "Grid generation successfully completed."

 stop


contains 

!*******************************************************************************
!* Prismatic Grid Generation
!*******************************************************************************
 subroutine prismatic_grid
 implicit none
!*******************************************************************************
! Generate prismatic grid
!*******************************************************************************
  ntrias_b = ntrias + ntrias        !Inner and outer boundaries
  allocate(tri(ntrias_b,5))

  if (n_sections==6) then
   nquads_b = nnodes_circum * (nr-1) !Symmetry boundary
   allocate(quad(nquads_b,5))
  endif

  nprs = ntrias*(nr-1)
  allocate(prs(nprs,6) )
  ntet = 0

  write(*,*) "1. Prismatic grid "
  write(*,*)
  write(*,*) "----- Predicted dimensions ------"
  write(*,*)
  write(*,*) "    nodes      = ", nnodes
  write(*,*) "  Volume elements:"
  write(*,*) "    prisms     = ", nprs
  write(*,*) "    tetrahedra = ", ntet
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b
  write(*,*) "         quads = ", nquads_b
  write(*,*)
  write(*,*) "---------------------------------"
  write(*,*)

  nprs = 0
  ntrias_b = 0
  nquads_b = 0

! Copy the triangulation on the body

  do i = 1, ntrias
          ntrias_b = ntrias_b + 1
   tri(ntrias_b,5) = 1
   tri(ntrias_b,1) = tria(i)%v(1)
   tri(ntrias_b,2) = tria(i)%v(2)
   tri(ntrias_b,3) = tria(i)%v(3)
  end do

! Generate prisms

  do i = 1, ntrias

   node1 = tria(i)%v(1)
   node2 = tria(i)%v(2)
   node3 = tria(i)%v(3)
   node4 = node_above(tria(i)%v(1))
   node5 = node_above(tria(i)%v(2))
   node6 = node_above(tria(i)%v(3))

     nprs = nprs + 1

     prs(nprs,1) = node1
     prs(nprs,2) = node2
     prs(nprs,3) = node3
     prs(nprs,4) = node4
     prs(nprs,5) = node5
     prs(nprs,6) = node6

   do k = 2, nr-1

    node1 = node4
    node2 = node5
    node3 = node6
    node4 = node_above(node1)
    node5 = node_above(node2)
    node6 = node_above(node3)

    nprs = nprs + 1
    prs(nprs,1) = node1
    prs(nprs,2) = node2
    prs(nprs,3) = node3
    prs(nprs,4) = node4
    prs(nprs,5) = node5
    prs(nprs,6) = node6


   end do

  !This is a triangle on the outer boundary

     ntrias_b = ntrias_b + 1
     tri(ntrias_b,1) = node6
     tri(ntrias_b,2) = node5
     tri(ntrias_b,3) = node4

     if ( node(node4)%x < -1.0e-13_dp .or. &
          node(node5)%x < -1.0e-13_dp .or. &
          node(node6)%x < -1.0e-13_dp      ) then
     !(1)Inflow
       tri(ntrias_b,5) = 2
     else
     !(2)Outflow
       tri(ntrias_b,5) = 3
     endif

  end do

  write(*,*) "------- Actual dimensions -------"
  write(*,*)
  write(*,*) "  Volume elements:"
  write(*,*) "    prisms     = ", nprs
  write(*,*) "    tetrahedra = ", ntet
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b
  write(*,*) "         quads = ", nquads_b
  write(*,*)
  write(*,*) "---------------------------------"
  write(*,*)

 end subroutine prismatic_grid


!*******************************************************************************
!* Tet grid Generation
!
!*******************************************************************************
 subroutine tet_grid
 implicit none
 integer :: tria_type
!*******************************************************************************
! Generate tet grid
!*******************************************************************************

  if (n_sections==6) then
   ntrias_b = ntrias + ntrias + 2*nnodes_circum * (nr-1) !Inner and outer boundaries
   allocate(tri(ntrias_b,5))
  elseif (n_sections==1) then
   ntrias_b = ntrias + ntrias + 2*nnodes_circum * (nr-1) !Inner and outer boundaries
   ntrias_b = ntrias_b + 2*2*(nr_gs)*(nr-1)
   allocate(tri(ntrias_b,5))
  endif

!  ntrias_b = ntrias + ntrias + 2*nnodes_circum * (nr-1) !Inner and outer boundaries
!  allocate(tri(ntrias_b,5))

  nquads_b = 0
  nprs = 0
  ntet = ntrias*(nr-1)*3
  allocate(tet(ntet,6) )

  write(*,*) "2. Tetrahedral grid "
  write(*,*)
  write(*,*) "----- Predicted dimensions ------"
  write(*,*)
  write(*,*) "    nodes      = ", nnodes
  write(*,*) "  Volume elements:"
  write(*,*) "    prisms     = ", nprs
  write(*,*) "    tetrahedra = ", ntet
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b
  write(*,*) "         quads = ", nquads_b
  write(*,*)
  write(*,*) "---------------------------------"
  write(*,*)

  nprs = 0
  ntet = 0
  ntrias_b = 0
  nquads_b = 0

! Copy the triangulation on the body: Reverse the orientation.

  do i = 1, ntrias
          ntrias_b = ntrias_b + 1
   tri(ntrias_b,5) = 1
   tri(ntrias_b,1) = tria(i)%v(1)
   tri(ntrias_b,2) = tria(i)%v(2)
   tri(ntrias_b,3) = tria(i)%v(3)
  end do

! Generate tetrahedra

  do i = 1, ntrias

   do k = 1, nr-1

    if (k == 1) then

     node1 = tria(i)%v(1)
     node2 = tria(i)%v(2)
     node3 = tria(i)%v(3)
     node4 = node_above(tria(i)%v(1))
     node5 = node_above(tria(i)%v(2))
     node6 = node_above(tria(i)%v(3))

     tria_type = tria(i)%type

    else

     node1 = node4
     node2 = node5
     node3 = node6
     node4 = node_above(node1)
     node5 = node_above(node2)
     node6 = node_above(node3)

    endif

    t_side : if (i < ntrias_hs+1) then

     t_type :if (tria_type == 1) then

      ntet = ntet + 1
      tet(ntet,1) = node1
      tet(ntet,2) = node2
      tet(ntet,3) = node3
      tet(ntet,4) = node4

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node6
      tet(ntet,3) = node2
      tet(ntet,4) = node4

      ntet = ntet + 1
      tet(ntet,1) = node6
      tet(ntet,2) = node3
      tet(ntet,3) = node2
      tet(ntet,4) = node4

     elseif (tria_type == 2) then

      ntet = ntet + 1
      tet(ntet,1) = node4
      tet(ntet,2) = node6
      tet(ntet,3) = node5
      tet(ntet,4) = node1

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node6
      tet(ntet,3) = node3
      tet(ntet,4) = node1

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node3
      tet(ntet,3) = node2
      tet(ntet,4) = node1

     endif t_type

    else

     t_type2 :if (tria_type == 1) then

      ntet = ntet + 1
      tet(ntet,1) = node1
      tet(ntet,2) = node2
      tet(ntet,3) = node3
      tet(ntet,4) = node4

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node6
      tet(ntet,3) = node3
      tet(ntet,4) = node4

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node3
      tet(ntet,3) = node2
      tet(ntet,4) = node4

     elseif (tria_type == 2) then

      ntet = ntet + 1
      tet(ntet,1) = node4
      tet(ntet,2) = node6
      tet(ntet,3) = node5
      tet(ntet,4) = node1

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node6
      tet(ntet,3) = node2
      tet(ntet,4) = node1

      ntet = ntet + 1
      tet(ntet,1) = node6
      tet(ntet,2) = node3
      tet(ntet,3) = node2
      tet(ntet,4) = node1

     endif t_type2

    endif t_side

   end do

  !This is a triangle on the outer boundary

     ntrias_b = ntrias_b + 1
     tri(ntrias_b,1) = node6
     tri(ntrias_b,2) = node5
     tri(ntrias_b,3) = node4

     if ( node(node4)%x < -1.0e-13_dp .or. &
          node(node5)%x < -1.0e-13_dp .or. &
          node(node6)%x < -1.0e-13_dp      ) then
     !(1)Inflow
       tri(ntrias_b,5) = 2
     else
     !(2)Outflow
       tri(ntrias_b,5) = 3
     endif

  end do

  write(*,*)
  write(*,*) "------- Actual dimensions -------"
  write(*,*)
  write(*,*) "    nodes      = ", nnodes
  write(*,*) "  Volume elements:"
  write(*,*) "    prisms     = ", nprs
  write(*,*) "    tetrahedra = ", ntet
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b
  write(*,*) "         quads = ", nquads_b
  write(*,*)
  write(*,*) "---------------------------------"
  write(*,*)

 end subroutine tet_grid

!*******************************************************************************
!* Mixed grid Generation
!
!*******************************************************************************
 subroutine mixed_grid
 implicit none
 integer :: tria_type
!*******************************************************************************
! Generate mixed grid
!*******************************************************************************

  ntrias_b = 2*ntrias + 2*nnodes_circum * (nr-1)-(nm-1) !All boundaries
!  allocate(tri(ntrias_b,5))
  nquads_b = nnodes_circum * (nm-1) !A part of symmetry boundary
!  allocate(quad(nquads_b,5))

  nprs = ntrias*(nm-1)              !Boundary layer
  allocate(prs(nprs,6) )
  ntet = ntrias*((nr-1)-(nm-1))*3   !Outer region
  allocate(tet(ntet,6) )

  allocate(tri(ntrias_b,5))
  allocate(quad(nquads_b,5))

  write(*,*) "3. Mixed grid "
  write(*,*)
  write(*,*) "----- Predicted dimensions ------"
  write(*,*)
  write(*,*) "    nodes      = ", nnodes
  write(*,*) "  Volume elements:"
  write(*,*) "    prisms     = ", nprs
  write(*,*) "    tetrahedra = ", ntet
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b
  write(*,*) "         quads = ", nquads_b
  write(*,*)
  write(*,*) "---------------------------------"
  write(*,*)

  nprs = 0
  ntet = 0
  ntrias_b = 0
  nquads_b = 0

! Copy the triangulation on the body

  do i = 1, ntrias
          ntrias_b = ntrias_b + 1
   tri(ntrias_b,5) = 1
   tri(ntrias_b,1) = tria(i)%v(1)
   tri(ntrias_b,2) = tria(i)%v(2)
   tri(ntrias_b,3) = tria(i)%v(3)
  end do

! Generate prisms

  do i = 1, ntrias

   node1 = tria(i)%v(1)
   node2 = tria(i)%v(2)
   node3 = tria(i)%v(3)
   node4 = node_above(tria(i)%v(1))
   node5 = node_above(tria(i)%v(2))
   node6 = node_above(tria(i)%v(3))
   tria_type = tria(i)%type

     nprs = nprs + 1

     prs(nprs,1) = node1
     prs(nprs,2) = node2
     prs(nprs,3) = node3
     prs(nprs,4) = node4
     prs(nprs,5) = node5
     prs(nprs,6) = node6

   do k = 2, nr-1

    node1 = node4
    node2 = node5
    node3 = node6
    node4 = node_above(node1)
    node5 = node_above(node2)
    node6 = node_above(node3)

    prs_or_tet : if (k < nm) then

     nprs = nprs + 1
     prs(nprs,1) = node1
     prs(nprs,2) = node2
     prs(nprs,3) = node3
     prs(nprs,4) = node4
     prs(nprs,5) = node5
     prs(nprs,6) = node6

    else


    t_side : if (i < ntrias_hs+1) then

     t_type :if (tria_type == 1) then

      ntet = ntet + 1
      tet(ntet,1) = node1
      tet(ntet,2) = node2
      tet(ntet,3) = node3
      tet(ntet,4) = node4

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node6
      tet(ntet,3) = node2
      tet(ntet,4) = node4

      ntet = ntet + 1
      tet(ntet,1) = node6
      tet(ntet,2) = node3
      tet(ntet,3) = node2
      tet(ntet,4) = node4

     elseif (tria_type == 2) then

      ntet = ntet + 1
      tet(ntet,1) = node4
      tet(ntet,2) = node6
      tet(ntet,3) = node5
      tet(ntet,4) = node1

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node6
      tet(ntet,3) = node3
      tet(ntet,4) = node1

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node3
      tet(ntet,3) = node2
      tet(ntet,4) = node1

     endif t_type

    else

     t_type2 :if (tria_type == 1) then

      ntet = ntet + 1
      tet(ntet,1) = node1
      tet(ntet,2) = node2
      tet(ntet,3) = node3
      tet(ntet,4) = node4

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node6
      tet(ntet,3) = node3
      tet(ntet,4) = node4

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node3
      tet(ntet,3) = node2
      tet(ntet,4) = node4

     elseif (tria_type == 2) then

      ntet = ntet + 1
      tet(ntet,1) = node4
      tet(ntet,2) = node6
      tet(ntet,3) = node5
      tet(ntet,4) = node1

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node6
      tet(ntet,3) = node2
      tet(ntet,4) = node1

      ntet = ntet + 1
      tet(ntet,1) = node6
      tet(ntet,2) = node3
      tet(ntet,3) = node2
      tet(ntet,4) = node1

     endif t_type2

    endif t_side

    endif prs_or_tet

   end do

  !This is a triangle on the outer boundary

     ntrias_b = ntrias_b + 1
     tri(ntrias_b,1) = node6
     tri(ntrias_b,2) = node5
     tri(ntrias_b,3) = node4

     if ( node(node4)%x < -1.0e-13_dp .or. &
          node(node5)%x < -1.0e-13_dp .or. &
          node(node6)%x < -1.0e-13_dp      ) then
     !(1)Inflow
       tri(ntrias_b,5) = 2
     else
     !(2)Outflow
       tri(ntrias_b,5) = 3
     endif

  end do

  write(*,*) "------- Actual dimensions -------"
  write(*,*)
  write(*,*) "  Volume elements:"
  write(*,*) "    prisms     = ", nprs
  write(*,*) "    tetrahedra = ", ntet
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b
  write(*,*) "         quads = ", nquads_b
  write(*,*)
  write(*,*) "---------------------------------"
  write(*,*)

 end subroutine mixed_grid

!*******************************************************************************
!* Write out a file containing line info (only within the viscous layer)):
!* This information is for the implicit line relaxations or line agglomeration.
!* This format is used by FUN3D.
!*******************************************************************************
 subroutine line_info_file(n_points)
 implicit none

 integer, intent(in) :: n_points

 integer :: node_below
 integer :: i, k, n_lines, n_total_points, max_points, min_points, i_count, inode,os

  open(unit=1, file=filename_lines, status="unknown", iostat=os)
  max_points = n_points
  min_points = n_points
  n_lines = nnodes_body
  n_total_points = n_lines * n_points

  write(*,*) " Total number of lines = ", n_lines
  write(*,*) "            max points = ", max_points
  write(*,*) "            min points = ", min_points

  write(1,*) n_lines, n_total_points, "Total lines and points"
  write(1,*) min_points, max_points, "Min and max points in line"

  i_count = 0

! Go up and write out the node info for every node on the body.
! NOTE: Go up to the given limit: n_points

  do i = 1, nnodes_body

     i_count = i_count + 1

!  1. First node (on the body)
     write(1,*) n_points, " Points in line for line = ", i_count
     write(1,'(i10,a10,3es30.20)') i, " x/y/z= ", node(i)%x, node(i)%y, node(i)%z
     node_below = i

!  2. Second node to the node below the last one.
   do k = 2, n_points-1
     write(1,'(i10)') node_above(node_below)
     node_below = node_above(node_below)
   end do

!  3. Last node
     inode = node_above(node_below)
     write(1,'(i10,a10,3es30.20)') inode, " x/y/z= ", node(inode)%x, node(inode)%y, node(inode)%z

  end do

 if (i_count /= n_lines) write(*,*) "Error: i_count /= n_lines"
 write(*,*)
 write(*,*) "lines_fmt file has been written: ", filename_lines

 close(1)
 end subroutine line_info_file

!*******************************************************************************
!* Write out a file containing line info (Lines go all the way up to the outer):
!* This information is for the implicit line relaxations or line agglomeration.
!* This format is used by FUN3D.
!*******************************************************************************
 subroutine line_info_file_all(n_points)
 implicit none

 integer, intent(in) :: n_points

 integer :: node_below
 integer :: i, k, n_lines, n_total_points, max_points, min_points, i_count, inode,os

  open(unit=1, file=filename_lines, status="unknown", iostat=os)
  max_points = n_points
  min_points = n_points
  n_lines = nnodes_body
  n_total_points = n_lines * n_points

  write(*,*) " Total number of lines = ", n_lines
  write(*,*) "            max points = ", max_points
  write(*,*) "            min points = ", min_points

  write(1,*) n_lines, n_total_points, "Total lines and points"
  write(1,*) min_points, max_points, "Min and max points in line"

  i_count = 0

! Go up and write out the node info for every node on the body.
! NOTE: Go up to the given limit: n_points

  do i = 1, nnodes_body

     i_count = i_count + 1

!  1. First node (on the body)
     write(1,*) n_points, " Points in line for line = ", i_count
     write(1,'(i10,a10,3es30.20)') i, " x/y/z= ", node(i)%x, node(i)%y, node(i)%z
     node_below = i

!  2. Second node to the node below the last one.
   do k = 2, n_points-1
     write(1,'(i10)') node_above(node_below)
     node_below = node_above(node_below)
   end do

!  3. Last node
     inode = node_above(node_below)
     write(1,'(i10,a10,3es30.20)') inode, " x/y/z= ", node(inode)%x, node(inode)%y, node(inode)%z

  end do

 if (i_count /= n_lines) write(*,*) "Error: i_count /= n_lines"
 write(*,*)
 write(*,*) "lines_fmt file has been written: ", filename_lines

 close(1)
 end subroutine line_info_file_all


!*******************************************************************************
! This subroutine writes a Tecplot file for the volume grid.
!*******************************************************************************
 subroutine write_tecplot_volume_file

 open(unit=8, file=filename_tecplot_v, status="unknown", iostat=os)
 write(8,*) 'TITLE = "GRID"'
 write(8,*) 'VARIABLES = "x","y","z","k1","k2","k3","k4"'

! Tetra Zone
  if (ntet > 0) then

   write(8,*) 'zone  n=', nnodes,',e=', ntet,' , et=tetrahedron, f=fepoint'
   do i = 1, nnodes
     write(8,'(3es20.10,4i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i)
   end do

   do i = 1, ntet
    write(8,'(4i10)') tet(i,1), tet(i,2), tet(i,3), tet(i,4)
   end do

  endif

! Prism zone
  if (nprs > 0) then

   write(8,*) 'zone  n=', nnodes,',e=', nprs,' , et=brick, f=fepoint'
   do i = 1, nnodes
     write(8,'(3es20.10,4i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i)
   end do

   do i = 1, nprs
    write(8,'(8i10)') prs(i,1), prs(i,2), prs(i,3), prs(i,3), &
                      prs(i,4), prs(i,5), prs(i,6), prs(i,6)
   end do

  endif

 close(8)

 end subroutine write_tecplot_volume_file

!*******************************************************************************
! This subroutine writes  a Tecplot file for boundaries.
!******************************************************************************
 subroutine write_tecplot_boundary_file

 integer :: ntrias_outer

 open(unit=7, file=filename_tecplot_b, status="unknown", iostat=os)
 write(7,*) 'TITLE = "GRID"'
 write(7,*) 'VARIABLES = "x","y","z","k1","k2","k3","k4"'

!-------------------------------------------------------
! Triangles on the sphere surface
!-------------------------------------------------------
 write(7,*) 'ZONE T="Body"  N=', nnodes,',E=', ntrias,' , ET=quadrilateral, F=FEPOINT'
  do i = 1, nnodes
    write(7,'(3ES20.10,4i13)') node(i)%x, node(i)%y, node(i)%z, &
                               k1(i),k2(i),k3(i),k4(i)
  end do
  do i = 1, ntrias
   write(7,'(4I10)') tri(i,1), tri(i,2), tri(i,3), tri(i,3)
  end do

!-------------------------------------------------------
! Triangles on the outer INFLOW boundary.
!-------------------------------------------------------
  ntrias_outer = 0
  do i = ntrias+1, ntrias_b
   if ( tri(i,5) == 2 ) ntrias_outer = ntrias_outer + 1
  end do

 write(7,*) 'ZONE T="Inflow"  N=', nnodes,',E=', ntrias_outer, &
            ' , ET=quadrilateral, F=FEPOINT'
  do i = 1, nnodes
    write(7,'(3ES20.10,4i13)') node(i)%x, node(i)%y, node(i)%z, &
                               k1(i),k2(i),k3(i),k4(i)
  end do
  do i = ntrias+1, ntrias_b
   if ( tri(i,5) == 2 ) write(7,'(4I10)') tri(i,1), tri(i,2), tri(i,3), tri(i,3)
  end do

!-------------------------------------------------------
! Triangles on the outer OUTFLOW boundary.
!-------------------------------------------------------
  ntrias_outer = 0
  do i = ntrias+1, ntrias_b
   if ( tri(i,5) == 3 ) ntrias_outer = ntrias_outer + 1
  end do

 write(7,*) 'ZONE T="Outflow"  N=', nnodes,',E=', ntrias_outer, &
            ' , ET=quadrilateral, F=FEPOINT'
  do i = 1, nnodes
    write(7,'(3ES20.10,4i13)') node(i)%x, node(i)%y, node(i)%z, &
                               k1(i),k2(i),k3(i),k4(i)
  end do
  do i = ntrias+1, ntrias_b
   if ( tri(i,5) == 3 ) write(7,'(4I10)') tri(i,1), tri(i,2), tri(i,3), tri(i,3)
  end do

 close(7)

 end subroutine write_tecplot_boundary_file

!*******************************************************************************
! This subroutine writes a ugrid file.
!*******************************************************************************
 subroutine write_ugrid_file

  if ( b8_ugrid_format ) then
    open(unit=9, file=filename_ugrid, form='unformatted',access="stream",&
                                      status='unknown', iostat=os )
    write(9) nnodes,   ntrias_b,    nquads_b,   ntet,    0, nprs, 0
  else
    open(unit=9, file=filename_ugrid, status="unknown", iostat=os)
    !                    #nodes, #tri_faces, #quad_faces, #tetra, #pyr, #prz,
    !                    #hex
    write(9,'(7I20)') nnodes,   ntrias_b,    nquads_b,   ntet,    0, nprs, 0
  endif

   if ( b8_ugrid_format ) then

  !FIXME: Gfortran returns an error at compilation (Line below)
  !   (( node(i)%x, node(i)%y, node(i)%z            ), i = 1, nnodes   ), &
  !               1
  ! Error: Expected a right parenthesis in expression at (1)
  !
  ! Let me know if you know the solution. Thanks.

   write(9)                                                            &
   (( node(i)%x, node(i)%y, node(i)%z            ), i = 1, nnodes   ), &
   (( tri(i,1), tri(i,2), tri(i,3)               ), i = 1, ntrias_b ), &
   (( quad(i,1), quad(i,2), quad(i,3), quad(i,4) ), i = 1, nquads_b ), &
   (( tri(i,5)                                   ), i = 1, ntrias_b ), &
   (( quad(i,5)                                  ), i = 1, nquads_b ), &
   (( tet(i,1), tet(i,2), tet(i,3), tet(i,4)     ), i = 1, ntet     ), &
   (( prs(i,1), prs(i,2), prs(i,3),                                    &
      prs(i,4), prs(i,5), prs(i,6)               ), i = 1, nprs     )

!  Code runs with the following. But FUN3D doesn't run with the resulting .ugrid.
!  write(9)                                                          &
!  ( node(i)%x, node(i)%y, node(i)%z            , i = 1, nnodes   ), &
!  ( tri(i,1), tri(i,2), tri(i,3)               , i = 1, ntrias_b ), &
!  ( quad(i,1), quad(i,2), quad(i,3), quad(i,4) , i = 1, nquads_b ), &
!  ( tri(i,5)                                   , i = 1, ntrias_b ), &
!  ( quad(i,5)                                  , i = 1, nquads_b ), &
!  ( tet(i,1), tet(i,2), tet(i,3), tet(i,4)     , i = 1, ntet     ), &
!  ( prs(i,1), prs(i,2), prs(i,3),                                   &
!    prs(i,4), prs(i,5), prs(i,6)               , i = 1, nprs     )
!
!   Below also doesn't work...
!    do i=1,nnodes
!     write(9) node(i)%x, node(i)%y, node(i)%z
!    end do
!    do i=1,ntrias_b
!     write(9) tri(i,1), tri(i,2), tri(i,3)
!    end do
!    do i=1,nquads_b
!     write(9) quad(i,1), quad(i,2), quad(i,3), quad(i,4)
!    end do
!    do i=1,ntrias_b
!     write(9) tri(i,5)
!    end do
!    do i=1,nquads_b
!     write(9) quad(i,5)
!    end do
!    do i=1,ntet
!     write(9) tet(i,1), tet(i,2), tet(i,3), tet(i,4)
!    end do
!    do i=1,nprs
!     write(9) prs(i,1), prs(i,2), prs(i,3),prs(i,4), prs(i,5), prs(i,6)
!    end do

   else

! Nodes
  do i = 1, nnodes
   write(9,'(3ES26.15)') node(i)%x, node(i)%y, node(i)%z
  end do

! Triangular faces = ntri
  if (ntrias_b > 0) then
   do i = 1, ntrias_b
    write(9,'(3I20)') tri(i,1), tri(i,2), tri(i,3)
   end do
  endif

! Quad faces = nquad
  if (nquads_b > 0) then
   do i = 1, nquads_b
    write(9,'(4I20)') quad(i,1), quad(i,2), quad(i,3), quad(i,4)
   end do
  endif

! Face tag
  if (ntrias_b > 0) then
   do i = 1, ntrias_b
    write(9,'(I110)')  tri(i,5)
   end do
  endif

  if (nquads_b > 0) then
   do i = 1, nquads_b
    write(9,'(I110)') quad(i,5)
   end do
  endif

! tet
  if (ntet > 0) then
   do i = 1, ntet
    write(9,'(4I20)') tet(i,1), tet(i,2), tet(i,3), tet(i,4)
   end do
  endif

! Prism
  if (nprs > 0) then
   do i = 1, nprs
    write(9,'(6I20)') prs(i,1), prs(i,2), prs(i,3), &
                      prs(i,4), prs(i,5), prs(i,6)
   end do
  endif

endif

  close(9)

 end subroutine write_ugrid_file


!*******************************************************************************
! This subroutine writes an k file
!******************************************************************************
 subroutine write_k_file

 open(unit=10, file=filename_k, status="unknown", iostat=os)

! Write node_number, k1, k2, k3, k4

  write(10,*) nnodes

 do i = 1, nnodes
  write(10,'(5i13)') i, k1(i),k2(i),k3(i),k4(i)
 end do

 close(10)

 end subroutine write_k_file


!********************************************************************************
!* This subroutine is useful to expand or shrink type(node_data_yz) arrays.
!*
!*  Array, x, will be allocated if the requested dimension is 1 (i.e., n=1)
!*                     expanded to the requested dimension, n, if n > dim(x).
!*                     shrunk to the requested dimension, n, if n < dim(x).
!*
!********************************************************************************
  subroutine my_alloc_ndxy_ptr(x,n)
  implicit none

  integer,                    intent(in   ) :: n
  type(node_data_xy), dimension(:), pointer :: x

  integer :: i
  type(node_data_xy), dimension(:), pointer :: temp

  if (n <= 0) then
   write(*,*) "my_alloc_ndxy_ptr received non-positive dimension. Stop."
   stop
  endif

! If initial, allocate and return
  if (.not.(associated(x))) then
   allocate(x(n))
   return
  endif

! If reallocation, create a pointer with a target of new dimension.
  allocate(temp(n))
    temp(n)%gnode = 0
    temp(n)%x     = zero
    temp(n)%y     = zero

! (1) Expand the array dimension
  if ( n > size(x) ) then

   do i = 1, size(x) 
    temp(i)%gnode = x(i)%gnode
    temp(i)%x     = x(i)%y
    temp(i)%y     = x(i)%z
  end do

! (2) Shrink the array dimension: the extra data, x(n+1:size(x)), discarded.
  else

   do i = 1, n
    temp(i)%gnode = x(i)%gnode
    temp(i)%y     = x(i)%y
    temp(i)%z     = x(i)%z
   end do

  endif

! Destroy the target of x
  deallocate(x)

! Re-assign the pointer
   x => temp

  return


  end subroutine my_alloc_ndxy_ptr
!********************************************************************************

function int_sign(n)
implicit none
integer, intent(in) :: n
integer :: int_sign

if (n >= 0) then
 int_sign = 1
else
 int_sign = -1
endif

end function int_sign



 function big_endian_io( opt_unit )

 integer, intent(in) :: opt_unit
 logical             :: big_endian_io
! one-byte integer
 integer, parameter :: i1 = selected_int_kind(2)
! two-byte integer
 integer, parameter :: i2 = selected_int_kind(4)
 integer(i1) :: byte_one, byte_two
! 00000000 00000001 big-endian binary
 integer(i2) :: two_byte_int = 1_i2

    open(opt_unit,status='scratch',form='unformatted')
      write( opt_unit) two_byte_int
      rewind(opt_unit)
      read(  opt_unit) byte_one, byte_two
    close(opt_unit)
    big_endian_io = ( byte_one == 0 .and. byte_two == 1 )

 end function big_endian_io



end program sphere_grid


