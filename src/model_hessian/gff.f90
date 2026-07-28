! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!! ========================================================================
!  This module implements various model Hessians, the actual definition
!  of the model Hessian is given in the description of each header, the
!  implementation following tries to do its very best to transform it back
!  from redundant internal coordinates (where the model is defined) to
!  Cartesian coordinates which are then used in any imaginable fashion
!  later on in the actual optimization.                        - SAW190131
!! ========================================================================

module xtb_modelhessian_gff
   implicit none

   public :: mh_gff
   private

contains

      subroutine mh_gff(Cart,nAtoms,Hess,at,s6,param,topo,neigh)
      use xtb_gfnff_data, only : TGFFData
      use xtb_gfnff_topology, only : TGFFTopology
      use xtb_gfnff_neighbor, only : TNeigh
      use xtb_type_timer
      use xtb_bmatrix, only : bmat_torsion
      use xtb_modelhessian_eeq, only : itabrow
      Implicit Integer (i-n)
      Implicit Real*8 (a-h, o-z)
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(in) :: topo
      type(TNeigh), intent(in) :: neigh

      Integer at(nAtoms)
      Real*8 Cart(3,nAtoms)
      Real*8 Hess((3*nAtoms)*(3*nAtoms+1)/2)
      Real*8 s6
!     charges qa taken from gffcom

      Real*8 rij(3),rjk(3),rkl(3),&
     &       si(3),sj(3),sk(3),&
     &       sl(3),sm(3),x(2),y(2),z(2),&
     &       xyz(3,4), C(3,4)
      Real*8 rAV(3,3), aAV(3,3),rkr, rkf
      Integer iANr(nAtoms)
      logical profile

      Data rAv/1.3500d+00,2.1000d+00,2.5300d+00,&
     &         2.1000d+00,2.8700d+00,3.4000d+00,&
     &         2.5300d+00,3.4000d+00,3.4000d+00/ 
      Data aAv/1.0000d+00,0.3949d+00,0.3949d+00,&
     &         0.3949d+00,0.2800d+00,0.2800d+00,&
     &         0.3949d+00,0.2800d+00,0.2800d+00/ 
!     Data rkr,rkf,rkt/0.4500D+00,0.1500D+00,0.5000D-02/ ! org
      Data rkr,rkf,rkt/0.4500D+00,0.3000D+00,0.75000  /  ! adjusted to account for redundant angles in old version

!cc VDWx-Parameters (Grimme) used for vdw-correction of model hessian 
      real*8 c6(55),c66
! NOT USED ANYMORE (BJ instead)
!     data vander &
! H, He
!    &     /0.91d0,0.92d0,&
! Li-Ne
!    &      0.75d0,1.28d0,1.35d0,1.32d0,1.27d0,1.22d0,1.17d0,1.13d0,&
! Na-Ar
!    &      1.04d0,1.24d0,1.49d0,1.56d0,1.55d0,1.53d0,1.49d0,1.45d0,&
! K, Ca
!    &      1.35d0,1.34d0,&
! Sc-Zn
!    &      1.42d0,1.42d0,1.42d0,1.42d0,1.42d0,&
!    &      1.42d0,1.42d0,1.42d0,1.42d0,1.42d0,&
! Ga-Kr
!    &      1.50d0,1.57d0,1.60d0,1.61d0,1.59d0,1.57d0,&
! Rb, Sr
!    &      1.48d0,1.46d0,&
! Y-Cd
!    &      1.49d0,1.49d0,1.49d0,1.49d0,1.49d0,&
!    &      1.49d0,1.49d0,1.49d0,1.49d0,1.49d0,&
! In, Sn, Sb, Te, I, Xe
!    &      1.52d0,1.64d0,1.71d0,1.72d0,1.72d0,1.71d0/ 

      Real*8 Zero, One, Two, Three, Four, Five, Six, Seven,&
     &       Eight, RNine, Ten, Half, Pi, SqrtP2, TwoP34,&
     &       TwoP54, One2C2
      Parameter(Zero =0.0D0, One  =1.0D0, Two=2.0D0, Three=3.0D0,&
     &          Four =4.0D0, Five =5.0D0, Six=6.0D0, Seven=7.0D0,&
     &          Eight=8.0D0, rNine=9.0D0, Ten=1.0D1, Half=0.5D0,&
     &          Pi    =3.141592653589793D0,&
     &          SqrtP2=0.8862269254527579D0,&
     &          TwoP34=0.2519794355383808D0,&
     &          TwoP54=5.914967172795612D0,&
     &          One2C2=0.2662567690426443D-04)

      type(tb_timer) :: timer
!     inline fct
      lina(i,j)=min(i,j)+max(i,j)*(max(i,j)-1)/2        
      ixyz(i,iAtom) = (iAtom-1)*3 + i
      Jnd(i,j) = i*(i-1)/2 +j
      Ind(i,iAtom,j,jAtom)=Jnd(Max(ixyz(i,iAtom),ixyz(j,jAtom)),&
     &                         Min(ixyz(i,iAtom),ixyz(j,jAtom)))

!     map heavy atoms to Z<=54
      do i=1,nAtoms
         iANr(i)=at(i)
         if(at(i).gt.54) iANr(i)=at(i) - 18
         if(at(i).gt.72) iANr(i)=at(i) - 32
         if(at(i).gt.56.and.at(i).lt.72) iANr(i)=39 ! set LNs to Y
         if(at(i).gt.86) iANr(i)=55 ! use same c6 for all atom types > 86
      enddo

      profile=.false.

      if (profile) call timer%new(3,.false.)
      if (profile) call timer%measure(1,'bond/vdw/qq')
!     D2 C6
      c6( 1)=0.14; c6( 2)=0.08; c6( 3)=1.61; c6( 4)=1.61; c6( 5)=3.13
      c6( 6)=1.75; c6( 7)=1.23; c6( 8)=0.70; c6( 9)=0.75; c6(10)=0.63
      c6(11)=5.71; c6(12)=5.71; c6(13)=10.79; c6(14)=9.23; c6(15)=7.84
      c6(16)=5.57; c6(17)=5.07; c6(18)=4.61; c6(19:30)=10.8; c6(31)=16.99
      c6(32)=17.10; c6(33)=16.37; c6(34)=12.64; c6(35)=12.47; c6(36)=12.01
      c6(37:48)=24.67; c6(49)=37.32; c6(50)=38.71; c6(51)=38.44
      c6(52)=31.74; c6(53)=31.50; c6(54)=29.99
      c6(55)=40.0

!hjw threshold reduced
      rZero=1.0d-10
      n3=3*nAtoms
     Hess = 0.0d0


!
!     Hessian for tension
!
      do ibond=1,neigh%nbond
         kAtom=neigh%blist(1,ibond)
         lAtom=neigh%blist(2,ibond)
         kr=itabrow(iANr(kAtom))
         lr=itabrow(iANr(lAtom))
            xkl=Cart(1,kAtom)-Cart(1,lAtom)
            ykl=Cart(2,kAtom)-Cart(2,lAtom)
            zkl=Cart(3,kAtom)-Cart(3,lAtom)
            rkl2 = xkl**2 + ykl**2 + zkl**2
            r0=rAv(kr,lr)
            alpha=aAv(kr,lr)
            gmm=Exp(-alpha*rkl2)*rkr*Exp(alpha*r0**2)
            Hxx=gmm*xkl*xkl/rkl2
            Hxy=gmm*xkl*ykl/rkl2
            Hxz=gmm*xkl*zkl/rkl2
            Hyy=gmm*ykl*ykl/rkl2
            Hyz=gmm*ykl*zkl/rkl2
            Hzz=gmm*zkl*zkl/rkl2
            Hess(Ind(1,kAtom,1,kAtom))=Hess(Ind(1,kAtom,1,kAtom))+Hxx
            Hess(Ind(2,kAtom,1,kAtom))=Hess(Ind(2,kAtom,1,kAtom))+Hxy
            Hess(Ind(2,kAtom,2,kAtom))=Hess(Ind(2,kAtom,2,kAtom))+Hyy
            Hess(Ind(3,kAtom,1,kAtom))=Hess(Ind(3,kAtom,1,kAtom))+Hxz
            Hess(Ind(3,kAtom,2,kAtom))=Hess(Ind(3,kAtom,2,kAtom))+Hyz
            Hess(Ind(3,kAtom,3,kAtom))=Hess(Ind(3,kAtom,3,kAtom))+Hzz
            Hess(Ind(1,kAtom,1,lAtom))=Hess(Ind(1,kAtom,1,lAtom))-Hxx
            Hess(Ind(1,kAtom,2,lAtom))=Hess(Ind(1,kAtom,2,lAtom))-Hxy
            Hess(Ind(1,kAtom,3,lAtom))=Hess(Ind(1,kAtom,3,lAtom))-Hxz
            Hess(Ind(2,kAtom,1,lAtom))=Hess(Ind(2,kAtom,1,lAtom))-Hxy
            Hess(Ind(2,kAtom,2,lAtom))=Hess(Ind(2,kAtom,2,lAtom))-Hyy
            Hess(Ind(2,kAtom,3,lAtom))=Hess(Ind(2,kAtom,3,lAtom))-Hyz
            Hess(Ind(3,kAtom,1,lAtom))=Hess(Ind(3,kAtom,1,lAtom))-Hxz
            Hess(Ind(3,kAtom,2,lAtom))=Hess(Ind(3,kAtom,2,lAtom))-Hyz
            Hess(Ind(3,kAtom,3,lAtom))=Hess(Ind(3,kAtom,3,lAtom))-Hzz
            Hess(Ind(1,lAtom,1,lAtom))=Hess(Ind(1,lAtom,1,lAtom))+Hxx
            Hess(Ind(2,lAtom,1,lAtom))=Hess(Ind(2,lAtom,1,lAtom))+Hxy
            Hess(Ind(2,lAtom,2,lAtom))=Hess(Ind(2,lAtom,2,lAtom))+Hyy
            Hess(Ind(3,lAtom,1,lAtom))=Hess(Ind(3,lAtom,1,lAtom))+Hxz
            Hess(Ind(3,lAtom,2,lAtom))=Hess(Ind(3,lAtom,2,lAtom))+Hyz
            Hess(Ind(3,lAtom,3,lAtom))=Hess(Ind(3,lAtom,3,lAtom))+Hzz
      End do
 
!
!     Hessian for Coulomb+dispersion
!

      do kAtom=1,nAtoms
         do lAtom=1,kAtom-1
            xkl=Cart(1,kAtom)-Cart(1,lAtom)
            ykl=Cart(2,kAtom)-Cart(2,lAtom)
            zkl=Cart(3,kAtom)-Cart(3,lAtom)
            rkl2 = xkl**2 + ykl**2 + zkl**2
            if(rkl2.gt.1600.d0) cycle
            c66=-s6*sqrt(c6(iANr(katom))*c6(iANr(latom))) ! D2 value
            rr=sqrt(rkl2)
            rr3=rr*rkl2
            r02=param%d3r0(lina(at(katom),at(latom)))
            rrpa=rr+sqrt(r02) ! qq damping with BJ radius
            cqq=2.0d0*topo%qa(kAtom)*topo%qa(lAtom) ! a bit upscaled
            call getqqxx(xkl,    cqq,c66,rr,rkl2,rr3,rrpa,r02,hxx)
            call getqqxy(xkl,ykl,cqq,c66,rr,rkl2,rr3,rrpa,r02,hxy)
            call getqqxy(xkl,zkl,cqq,c66,rr,rkl2,rr3,rrpa,r02,hxz)
            call getqqxx(ykl,    cqq,c66,rr,rkl2,rr3,rrpa,r02,hyy)
            call getqqxy(ykl,zkl,cqq,c66,rr,rkl2,rr3,rrpa,r02,hyz)
            call getqqxx(zkl,    cqq,c66,rr,rkl2,rr3,rrpa,r02,hzz)
            Hess(Ind(1,kAtom,1,kAtom))=Hess(Ind(1,kAtom,1,kAtom))+Hxx
            Hess(Ind(2,kAtom,1,kAtom))=Hess(Ind(2,kAtom,1,kAtom))+Hxy
            Hess(Ind(2,kAtom,2,kAtom))=Hess(Ind(2,kAtom,2,kAtom))+Hyy
            Hess(Ind(3,kAtom,1,kAtom))=Hess(Ind(3,kAtom,1,kAtom))+Hxz
            Hess(Ind(3,kAtom,2,kAtom))=Hess(Ind(3,kAtom,2,kAtom))+Hyz
            Hess(Ind(3,kAtom,3,kAtom))=Hess(Ind(3,kAtom,3,kAtom))+Hzz
            Hess(Ind(1,kAtom,1,lAtom))=Hess(Ind(1,kAtom,1,lAtom))-Hxx
            Hess(Ind(1,kAtom,2,lAtom))=Hess(Ind(1,kAtom,2,lAtom))-Hxy
            Hess(Ind(1,kAtom,3,lAtom))=Hess(Ind(1,kAtom,3,lAtom))-Hxz
            Hess(Ind(2,kAtom,1,lAtom))=Hess(Ind(2,kAtom,1,lAtom))-Hxy
            Hess(Ind(2,kAtom,2,lAtom))=Hess(Ind(2,kAtom,2,lAtom))-Hyy
            Hess(Ind(2,kAtom,3,lAtom))=Hess(Ind(2,kAtom,3,lAtom))-Hyz
            Hess(Ind(3,kAtom,1,lAtom))=Hess(Ind(3,kAtom,1,lAtom))-Hxz
            Hess(Ind(3,kAtom,2,lAtom))=Hess(Ind(3,kAtom,2,lAtom))-Hyz
            Hess(Ind(3,kAtom,3,lAtom))=Hess(Ind(3,kAtom,3,lAtom))-Hzz
            Hess(Ind(1,lAtom,1,lAtom))=Hess(Ind(1,lAtom,1,lAtom))+Hxx
            Hess(Ind(2,lAtom,1,lAtom))=Hess(Ind(2,lAtom,1,lAtom))+Hxy
            Hess(Ind(2,lAtom,2,lAtom))=Hess(Ind(2,lAtom,2,lAtom))+Hyy
            Hess(Ind(3,lAtom,1,lAtom))=Hess(Ind(3,lAtom,1,lAtom))+Hxz
            Hess(Ind(3,lAtom,2,lAtom))=Hess(Ind(3,lAtom,2,lAtom))+Hyz
            Hess(Ind(3,lAtom,3,lAtom))=Hess(Ind(3,lAtom,3,lAtom))+Hzz
         End Do
      End Do

      if (profile) call timer%measure(1)

!
!     Hessian for bending
!
      if (profile) call timer%measure(2,'bend')
      do iangl=1,topo%nangl
         mAtom=topo%alist(1,iangl)
         iAtom=topo%alist(2,iangl)
         jAtom=topo%alist(3,iangl)
         mr=itabrow(iANr(mAtom))
         ir=itabrow(iANr(iAtom))
         jr=itabrow(iANr(jAtom))
 
            xmi=(Cart(1,iAtom)-Cart(1,mAtom))
            ymi=(Cart(2,iAtom)-Cart(2,mAtom))
            zmi=(Cart(3,iAtom)-Cart(3,mAtom))
            rmi2 = xmi**2 + ymi**2 + zmi**2
            rmi=sqrt(rmi2)
            r0mi=rAv(mr,ir)
            ami=aAv(mr,ir)
            xmj=(Cart(1,jAtom)-Cart(1,mAtom))
            ymj=(Cart(2,jAtom)-Cart(2,mAtom))
            zmj=(Cart(3,jAtom)-Cart(3,mAtom))
            rmj2 = xmj**2 + ymj**2 + zmj**2
            rmj=sqrt(rmj2)
            r0mj=rAv(mr,jr)
            amj=aAv(mr,jr)
            xij=(Cart(1,jAtom)-Cart(1,iAtom))
            yij=(Cart(2,jAtom)-Cart(2,iAtom))
            zij=(Cart(3,jAtom)-Cart(3,iAtom))
            rij2 = xij**2 + yij**2 + zij**2
            rrij=sqrt(rij2)

            alpha=rkf*exp((ami*r0mi**2+amj*r0mj**2))
            r=sqrt(rmj2+rmi2)
            gij=alpha*exp(-(ami*rmi2+amj*rmj2))
            rL2=(ymi*zmj-zmi*ymj)**2+(zmi*xmj-xmi*zmj)**2+&
     &         (xmi*ymj-ymi*xmj)**2
!hjw modified
            if(rL2.lt.1.d-14) then
              rL=0
            else
              rL=sqrt(rL2)
            end if
            if ((rmj.gt.rZero).and.(rmi.gt.rZero).and.&
     &                                (rrij.gt.rZero)) Then
              SinPhi=rL/(rmj*rmi)
              rmidotrmj=xmi*xmj+ymi*ymj+zmi*zmj
              CosPhi=rmidotrmj/(rmj*rmi)
!-------------None linear case
              If (SinPhi.gt.rZero) Then
                si(1)=(xmi/rmi*cosphi-xmj/rmj)/(rmi*sinphi)
                si(2)=(ymi/rmi*cosphi-ymj/rmj)/(rmi*sinphi)
                si(3)=(zmi/rmi*cosphi-zmj/rmj)/(rmi*sinphi)
                sj(1)=(cosphi*xmj/rmj-xmi/rmi)/(rmj*sinphi)
                sj(2)=(cosphi*ymj/rmj-ymi/rmi)/(rmj*sinphi)
                sj(3)=(cosphi*zmj/rmj-zmi/rmi)/(rmj*sinphi)
                sm(1)=-si(1)-sj(1)
                sm(2)=-si(2)-sj(2)
                sm(3)=-si(3)-sj(3)
                Do icoor=1,3
                   Do jCoor=1,3
                    If (mAtom.gt.iAtom) Then
                       Hess(Ind(icoor,mAtom,jcoor,iAtom))=&
     &                        Hess(Ind(icoor,mAtom,jcoor,iAtom))&
     &                        +gij*sm(icoor)*si(jcoor)
                    else
                      Hess(Ind(icoor,iAtom,jcoor,mAtom))=&
     &                        Hess(Ind(icoor,iAtom,jcoor,mAtom))&
     &                        +gij*si(icoor)*sm(jcoor)
                    End If
                    If (mAtom.gt.jAtom) Then
                        Hess(Ind(icoor,mAtom,jcoor,jAtom))=&
     &                        Hess(Ind(icoor,mAtom,jcoor,jAtom))&
     &                        +gij*sm(icoor)*sj(jcoor)
                    else
                      Hess(Ind(icoor,jAtom,jcoor,mAtom))=&
     &                        Hess(Ind(icoor,jAtom,jcoor,mAtom))&
     &                        +gij*sj(icoor)*sm(jcoor)
                    End If
                    If (iAtom.gt.jAtom) Then
                        Hess(Ind(icoor,iAtom,jcoor,jAtom))=&
     &                        Hess(Ind(icoor,iAtom,jcoor,jAtom))&
     &                        +gij*si(icoor)*sj(jcoor)
                     else
                        Hess(Ind(icoor,jAtom,jcoor,iAtom))=&
     &                        Hess(Ind(icoor,jAtom,jcoor,iAtom))&
     &                        +gij*sj(icoor)*si(jcoor)
                     End If
                   End Do
                End Do
                Do icoor=1,3
                  Do jCoor=1,icoor
                    Hess(Ind(icoor,iAtom,jcoor,iAtom))=&
     &                        Hess(Ind(icoor,iAtom,jcoor,iAtom))&
     &                        +gij*si(icoor)*si(jcoor)
                    Hess(Ind(icoor,mAtom,jcoor,mAtom))=&
     &                        Hess(Ind(icoor,mAtom,jcoor,mAtom))&
     &                        +gij*sm(icoor)*sm(jcoor)
                    Hess(Ind(icoor,jAtom,jcoor,jAtom))=&
     &                        Hess(Ind(icoor,jAtom,jcoor,jAtom))&
     &                        +gij*sj(icoor)*sj(jcoor)

                  End Do
                End Do
              Else
!----------------Linear case
                    if ((abs(ymi).gt.rZero).or.(abs(xmi).gt.rZero)) Then
                      x(1)=-ymi
                      y(1)=xmi
                      z(1)=Zero
                      x(2)=-xmi*zmi
                      y(2)=-ymi*zmi
                      z(2)=xmi*xmi+ymi*ymi
                    Else
                      x(1)=One
                      y(1)=Zero
                      z(1)=Zero
                      x(2)=Zero
                      y(2)=One
                      z(2)=Zero
                    End If
                    Do i=1,2
                     r1=sqrt(x(i)**2+y(i)**2+z(i)**2)
                     cosThetax=x(i)/r1
                     cosThetay=y(i)/r1
                     cosThetaz=z(i)/r1
                     si(1)=-cosThetax/rmi
                     si(2)=-cosThetay/rmi
                     si(3)=-cosThetaz/rmi
                     sj(1)=-cosThetax/rmj
                     sj(2)=-cosThetay/rmj
                     sj(3)=-cosThetaz/rmj
                     sm(1)=-(si(1)+sj(1))
                     sm(2)=-(si(2)+sj(2))
                     sm(3)=-(si(3)+sj(3))
                     Do icoor=1,3
                       Do jCoor=1,3
                        If (mAtom.gt.iAtom) Then
                          Hess(Ind(icoor,mAtom,jcoor,iAtom))=&
     &                        Hess(Ind(icoor,mAtom,jcoor,iAtom))&
     &                         +gij*sm(icoor)*si(jcoor)
                        else
                           Hess(Ind(icoor,iAtom,jcoor,mAtom))=&
     &                        Hess(Ind(icoor,iAtom,jcoor,mAtom))&
     &                         +gij*si(icoor)*sm(jcoor)
                        End If
                        If (mAtom.gt.jAtom) Then
                          Hess(Ind(icoor,mAtom,jcoor,jAtom))=&
     &                        Hess(Ind(icoor,mAtom,jcoor,jAtom))&
     &                         +gij*sm(icoor)*sj(jcoor)
                        else
                          Hess(Ind(icoor,jAtom,jcoor,mAtom))=&
     &                        Hess(Ind(icoor,jAtom,jcoor,mAtom))&
     &                         +gij*sj(icoor)*sm(jcoor)
                        End If
                        If (iAtom.gt.jAtom) Then
                           Hess(Ind(icoor,iAtom,jcoor,jAtom))=&
     &                        Hess(Ind(icoor,iAtom,jcoor,jAtom))&
     &                         +gij*si(icoor)*sj(jcoor)
                        else
                           Hess(Ind(icoor,jAtom,jcoor,iAtom))=&
     &                        Hess(Ind(icoor,jAtom,jcoor,iAtom))&
     &                         +gij*sj(icoor)*si(jcoor)
                        End If
                       End Do
                     End Do
                     Do icoor=1,3
                       Do jCoor=1,icoor
                         Hess(Ind(icoor,iAtom,jcoor,iAtom))=&
     &                        Hess(Ind(icoor,iAtom,jcoor,iAtom))&
     &                         +gij*si(icoor)*si(jcoor)
                         Hess(Ind(icoor,mAtom,jcoor,mAtom))=&
     &                        Hess(Ind(icoor,mAtom,jcoor,mAtom))&
     &                         +gij*sm(icoor)*sm(jcoor)
                         Hess(Ind(icoor,jAtom,jcoor,jAtom))=&
     &                         Hess(Ind(icoor,jAtom,jcoor,jAtom))&
     &                         +gij*sj(icoor)*sj(jcoor)
                       End Do
                     End Do
                 End Do
              End If
            End If
      End Do ! bend list

      if (profile) call timer%measure(2)

!
!     Hessian for torsion
!
      if (profile) call timer%measure(3,'torsion')
      do itors=1,topo%ntors
         iAtom=topo%tlist(3,itors)
         jAtom=topo%tlist(1,itors)
         kAtom=topo%tlist(2,itors)
         lAtom=topo%tlist(4,itors)

                 ir=itabrow(iANr(iAtom))
                 jr=itabrow(iANr(jAtom))
                 kr=itabrow(iANr(kAtom))
                 lr=itabrow(iANr(lAtom))
                 xyz(1:3,1)=Cart(1:3,iAtom)
                 xyz(1:3,2)=Cart(1:3,jAtom)
                 xyz(1:3,3)=Cart(1:3,kAtom)
                 xyz(1:3,4)=Cart(1:3,lAtom)
                 rij(1)=Cart(1,iAtom)-Cart(1,jAtom)
                 rij(2)=Cart(2,iAtom)-Cart(2,jAtom)
                 rij(3)=Cart(3,iAtom)-Cart(3,jAtom)
                 rij0=rAv(ir,jr)**2
                 aij =aAv(ir,jr)
                 rjk(1)=Cart(1,jAtom)-Cart(1,kAtom)
                 rjk(2)=Cart(2,jAtom)-Cart(2,kAtom)
                 rjk(3)=Cart(3,jAtom)-Cart(3,kAtom)
                 rjk0=rAv(jr,kr)**2
                 ajk =aAv(jr,kr)
                 rkl(1)=Cart(1,kAtom)-Cart(1,lAtom)
                 rkl(2)=Cart(2,kAtom)-Cart(2,lAtom)
                 rkl(3)=Cart(3,kAtom)-Cart(3,lAtom)
                 rkl0=rAv(kr,lr)**2
                 akl =aAv(kr,lr)
                 rij2=rij(1)**2+rij(2)**2+rij(3)**2
                 rjk2=rjk(1)**2+rjk(2)**2+rjk(3)**2
                 rkl2=rkl(1)**2+rkl(2)**2+rkl(3)**2

               beta=rkt*exp( (aij*rij0+ajk*rjk0+akl*rkl0))
               tij=beta*exp(-(aij*rij2+ajk*rjk2+akl*rkl2))

               C = bmat_torsion(xyz)
               si(1:3)=C(1:3,1)
               sj(1:3)=C(1:3,2)
               sk(1:3)=C(1:3,3)
               sl(1:3)=C(1:3,4)
!-------------Off diagonal block
              Do icoor=1,3
                Do jCoor=1,3
                 Hess(Ind(icoor,iAtom,jcoor,jAtom))=&
     &           Hess(Ind(icoor,iAtom,jcoor,jAtom))&
     &            +tij*si(icoor) * sj(jcoor)
                 Hess(Ind(icoor,iAtom,jcoor,kAtom))=&
     &           Hess(Ind(icoor,iAtom,jcoor,kAtom))&
     &            +tij*si(icoor) * sk(jcoor)
                 Hess(Ind(icoor,iAtom,jcoor,lAtom))=&
     &           Hess(Ind(icoor,iAtom,jcoor,lAtom))&
     &            +tij*si(icoor) * sl(jcoor)
                 Hess(Ind(icoor,jAtom,jcoor,kAtom))=&
     &           Hess(Ind(icoor,jAtom,jcoor,kAtom))&
     &            +tij*sj(icoor) * sk(jcoor)
                 Hess(Ind(icoor,jAtom,jcoor,lAtom))=&
     &           Hess(Ind(icoor,jAtom,jcoor,lAtom))&
     &            +tij*sj(icoor) * sl(jcoor)
                 Hess(Ind(icoor,kAtom,jcoor,lAtom))=&
     &           Hess(Ind(icoor,kAtom,jcoor,lAtom))&
     &            +tij*sk(icoor) * sl(jcoor)

                End Do
              End Do
!-------------Diagonal block
              Do icoor=1,3
                Do jCoor=1,icoor
                 Hess(Ind(icoor,iAtom,jcoor,iAtom))=&
     &           Hess(Ind(icoor,iAtom,jcoor,iAtom))&
     &            +tij*si(icoor) * si(jcoor)
                 Hess(Ind(icoor,jAtom,jcoor,jAtom))=&
     &           Hess(Ind(icoor,jAtom,jcoor,jAtom))&
     &            +tij*sj(icoor) * sj(jcoor)
                 Hess(Ind(icoor,kAtom,jcoor,kAtom))=&
     &           Hess(Ind(icoor,kAtom,jcoor,kAtom))&
     &            +tij*sk(icoor) * sk(jcoor)
                 Hess(Ind(icoor,lAtom,jcoor,lAtom))=&
     &           Hess(Ind(icoor,lAtom,jcoor,lAtom))&
     &            +tij*sl(icoor) * sl(jcoor)

                 End Do
               End Do
 
      enddo ! tors list

      if (profile) call timer%measure(3)
      if (profile) call timer%write(6,'modhes')
      Return
      end subroutine mh_gff

! Coulomb + dispersion Hessian (SG, 5/19)
! d/(dxdx)  q*q/(damp+r),  rpa=damp+r
! d/(dxdx)   c6/(r02^3+r^6)
      subroutine getqqxx(dx,qq,c6,r,r2,r3,rpa,r02,d2)      
      implicit none
      real*8 dx,qq,c6,r,r2,r3,rpa,r02,d2
      real*8 rpa2,dx2,ar6,r6,r8

!     Coulomb part
      rpa2=rpa**2
      dx2 = dx**2
      d2  = qq*(2.d0*dx2/(r2*rpa*rpa2) + dx2/(r3*rpa2) - 1.d0/(r*rpa2))
!     disp part
      r6  = r3*r3
      r8  = r6*r2
      ar6 = r02**3+r6
      d2  = d2 + c6*(dx2*72.d0*r8/ar6**3 - dx2*24.d0*r2/ar6**2 &
     &         - 6.d0*r2*r2/ar6**2)
      end subroutine getqqxx

! d/(dxdy)  q*q/(damp+r),  rpa=damp+r
! d/(dxdy)   c6/(r02^3+r^6)
      subroutine getqqxy(dx,dy,qq,c6,r,r2,r3,rpa,r02,d2)      
      implicit none
      real*8 dx,dy,qq,c6,r,r2,r3,rpa,r02,d2
      real*8 rpa2,r6,r8,ar6

!     Coulomb part
      rpa2=rpa**2
      d2  =qq*(2.d0*dx*dy/(r2*rpa*rpa2) + dx*dy/(r3*rpa2))
!     disp part
      r6  = r3*r3
      r8  = r6*r2
      ar6 = r02**3+r6
      d2  = d2 + c6*(dx*dy*72.d0*r8/ar6**3 - dx*dy*24.d0*r2/ar6**2)

      end subroutine getqqxy

end module xtb_modelhessian_gff
