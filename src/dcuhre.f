c
c     this is a summary of the routines in this file. 
c
c     dtest1   simple test driver for dcuhre.
c              sample output from a sun 3/50 is included.
c     dcuhre   main integrator. dcuhre calls dchhre and dadhre.
c     dchhre   checks the input to dcuhre.
c     dadhre   the adaptive integration routine. 
c              dadhre calls dtrhre, dinhre and drlhre.
c     dtrhre   maintaines the heap of subregions.
c     dinhre   computes weights and abscissas of the integration
c              rule. dinhre calls d132re, d112re, d09hre and d07hre.
c     d132re   computes weights and abscissas for a 2-dimensional
c              rule of degree 13.
c     d113re   computes weights and abscissas for a 3-dimensional 
c              rule of degree 11.
c     d09hre   computes weights and abscissas for a degree 9 rule.
c     d07hre   computes weights and abscissas for a degree 7 rule.
c     drlhre   computes estimates of integral and error over
c              subregions.
c     dfshre   computes fully symmetric sums of function evaluations.
c
c   dtest1 is a simple test driver for dcuhre.
c
c   output produced on a sun 3/50.
c
c       dcuhre test results
c
c    ftest calls = 3549, ifail =  0
c   n   estimated error    integral
c   1       0.000000       0.138508
c   2       0.000000       0.063695
c   3       0.000009       0.058618
c   4       0.000000       0.054070
c   5       0.000000       0.050056
c   6       0.000000       0.046546
c
c$$$      program dtest1
c$$$      external ftest
c$$$      integer key, n, nf, ndim, mincls, maxcls, ifail, neval, nw
c$$$      parameter (ndim = 4, nw = 5000, nf = 1)
c$$$      double precision a(ndim), b(ndim), wrkstr(nw)
c$$$      double precision absest(nf), finest(nf), absreq, relreq
c$$$      do 10 n = 1,ndim
c$$$         a(n) = 0
c$$$         b(n) = 3.141592
c$$$   10 continue
c$$$      mincls = 0
c$$$      maxcls = 50000
c$$$      key    = 0
c$$$      absreq = 0
c$$$      relreq = 1e-5
c$$$      call dcuhre(ndim, nf, a, b, mincls, maxcls, ftest, absreq, relreq,
c$$$     * key, nw, 0, finest, absest, neval, ifail, wrkstr)
c$$$      print 9999, neval, ifail
c$$$ 9999 format (8x, 'dcuhre test results', //'     ftest calls = ', i6,
c$$$     * ', ifail = ', i2, /'    n   estimated error    integral')
c$$$      do 20 n = 1,nf
c$$$         print 9998, n, absest(n), finest(n)
c$$$ 9998    format (3x, i2, 2f15.6)
c$$$   20 continue
c$$$      end
c$$$
c$$$c==============================================================
c$$$
c$$$
c$$$      subroutine ftest(ndim, z, nfun, f)
c$$$      integer n, ndim, nfun
c$$$      double precision z(ndim), f(nfun), sum
c$$$      sum = 0
c$$$c      do 10 n = 1,ndim
c$$$c         sum = sum + n*z(n)**2
c$$$c   10 continue
c$$$c      f(1) = exp(-sum/2)
c$$$c      f(1) = sin(z(1)*z(2))+sin(z(3)*z(4))
c$$$c      f(1) = z(1)*z(2)*z(3)*z(4)
c$$$       call f0(ndim, z, f(1))
c$$$c      do 20 n = 1,ndim
c$$$c         f(1) = z(n)*f(1)
c$$$c   20 continue
c$$$      end
c$$$
c$$$c=================================================================
c$$$      
c$$$      subroutine f0(ndim, z, f)
c$$$      integer          ndim
c$$$      double precision z(ndim), f
c$$$
c$$$      f = sin(z(1)*z(2))+sin(z(3)*z(4))
c$$$      
c$$$      end
c$$$

c=============================================================================

      subroutine dcuhre(ndim,numfun,a,b,minpts,maxpts,funsub,epsabs,
     +                  epsrel,key,nw,restar,result,abserr,neval,ifail,
     +                  work)
c***begin prologue dcuhre
c***date written   900116   (yymmdd)
c***revision date  900116   (yymmdd)
c***category no. h2b1a1
c***author
c            jarle berntsen, the computing centre,
c            university of bergen, thormohlens gt. 55,
c            n-5008 bergen, norway
c            phone..  47-5-544055
c            email..  jarle@eik.ii.uib.no
c            terje o. espelid, department of informatics,
c            university of bergen, thormohlens gt. 55,
c            n-5008 bergen, norway
c            phone..  47-5-544180
c            email..  terje@eik.ii.uib.no
c            alan genz, computer science department, washington state
c            university, pullman, wa 99163-2752, usa
c            email..  acg@eecs.wsu.edu
c***keywords automatic multidimensional integrator,
c            n-dimensional hyper-rectangles,
c            general purpose, global adaptive
c***purpose  the routine calculates an approximation to a given
c            vector of definite integrals
c
c      b(1) b(2)     b(ndim)
c     i    i    ... i       (f ,f ,...,f      ) dx(ndim)...dx(2)dx(1),
c      a(1) a(2)     a(ndim)  1  2      numfun
c
c       where f = f (x ,x ,...,x    ), i = 1,2,...,numfun.
c              i   i  1  2      ndim
c
c            hopefully satisfying for each component of i the following
c            claim for accuracy:
c            abs(i(k)-result(k)).le.max(epsabs,epsrel*abs(i(k)))
c***description computation of integrals over hyper-rectangular
c            regions.
c            dcuhre is a driver for the integration routine
c            dadhre, which repeatedly subdivides the region
c            of integration and estimates the integrals and the
c            errors over the subregions with greatest
c            estimated errors until the error request
c            is met or maxpts function evaluations have been used.
c
c            for ndim = 2 the default integration rule is of 
c            degree 13 and uses 65 evaluation points.
c            for ndim = 3 the default integration rule is of 
c            degree 11 and uses 127 evaluation points.
c            for ndim greater then 3 the default integration rule
c            is of degree 9 and uses num evaluation points where
c              num = 1 + 4*2*ndim + 2*ndim*(ndim-1) + 4*ndim*(ndim-1) +
c                    4*ndim*(ndim-1)*(ndim-2)/3 + 2**ndim
c            the degree 9 rule may also be applied for ndim = 2
c            and ndim = 3.
c            a rule of degree 7 is available in all dimensions.
c            the number of evaluation
c            points used by the degree 7 rule is
c              num = 1 + 3*2*ndim + 2*ndim*(ndim-1) + 2**ndim
c
c            when dcuhre computes estimates to a vector of
c            integrals, all components of the vector are given
c            the same treatment. that is, i(f ) and i(f ) for
c                                            j         k 
c            j not equal to k, are estimated with the same
c            subdivision of the region of integration.
c            for integrals with enough similarity, we may save
c            time by applying dcuhre to all integrands in one call.
c            for integrals that varies continuously as functions of
c            some parameter, the estimates produced by dcuhre will
c            also vary continuously when the same subdivision is
c            applied to all components. this will generally not be
c            the case when the different components are given
c            separate treatment.
c
c            on the other hand this feature should be used with
c            caution when the different components of the integrals
c            require clearly different subdivisions.
c
c   on entry
c
c     ndim   integer.
c            number of variables. 1 < ndim <=  15.
c     numfun integer.
c            number of components of the integral.
c     a      real array of dimension ndim.
c            lower limits of integration.
c     b      real array of dimension ndim.
c            upper limits of integration.
c     minpts integer.
c            minimum number of function evaluations.
c     maxpts integer.
c            maximum number of function evaluations.
c            the number of function evaluations over each subregion
c            is num.
c            if (key = 0 or key = 1) and (ndim = 2) then
c              num = 65
c            elseif (key = 0 or key = 2) and (ndim = 3) then
c              num = 127
c            elseif (key = 0 and ndim > 3) or (key = 3) then
c              num = 1 + 4*2*ndim + 2*ndim*(ndim-1) + 4*ndim*(ndim-1) +
c                    4*ndim*(ndim-1)*(ndim-2)/3 + 2**ndim
c            elseif (key = 4) then
c              num = 1 + 3*2*ndim + 2*ndim*(ndim-1) + 2**ndim
c            maxpts >= 3*num and maxpts >= minpts
c            for 3 < ndim < 13 the minimum values for maxpts are:
c             ndim =    4   5   6    7    8    9    10   11    12
c            key = 3:  459 819 1359 2151 3315 5067 7815 12351 20235
c            key = 4:  195 309  483  765 1251 2133 3795  7005 13299
c     funsub externally declared subroutine for computing
c            all components of the integrand at the given
c            evaluation point.
c            it must have parameters (ndim,x,numfun,funvls)
c            input parameters:
c              ndim   integer that defines the dimension of the
c                     integral.
c              x      real array of dimension ndim
c                     that defines the evaluation point.
c              numfun integer that defines the number of
c                     components of i.
c            output parameter:
c              funvls real array of dimension numfun
c                     that defines numfun components of the integrand.
c
c     epsabs real.
c            requested absolute error.
c     epsrel real.
c            requested relative error.
c     key    integer.
c            key to selected local integration rule.
c            key = 0 is the default value.
c                  for ndim = 2 the degree 13 rule is selected.
c                  for ndim = 3 the degree 11 rule is selected.
c                  for ndim > 3 the degree  9 rule is selected.
c            key = 1 gives the user the 2 dimensional degree 13
c                  integration rule that uses 65 evaluation points.
c            key = 2 gives the user the 3 dimensional degree 11
c                  integration rule that uses 127 evaluation points.
c            key = 3 gives the user the degree 9 integration rule.
c            key = 4 gives the user the degree 7 integration rule.
c                  this is the recommended rule for problems that
c                  require great adaptivity.
c     nw     integer.
c            defines the length of the working array work.
c            let maxsub denote the maximum allowed number of subregions
c            for the given values of maxpts, key and ndim.
c            maxsub = (maxpts-num)/(2*num) + 1
c            nw should be greater or equal to
c            maxsub*(2*ndim+2*numfun+2) + 17*numfun + 1
c            for efficient execution on parallel computers
c            nw should be chosen greater or equal to
c            maxsub*(2*ndim+2*numfun+2) + 17*numfun*mdiv + 1
c            where mdiv is the number of subregions that are divided in
c            each subdivision step.
c            mdiv is default set internally in dcuhre equal to 1.
c            for efficient execution on parallel computers
c            with nproc processors mdiv should be set equal to
c            the smallest integer such that mod(2*mdiv,nproc) = 0.
c
c     restar integer.
c            if restar = 0, this is the first attempt to compute
c            the integral.
c            if restar = 1, then we restart a previous attempt.
c            in this case the only parameters for dcuhre that may
c            be changed (with respect to the previous call of dcuhre)
c            are minpts, maxpts, epsabs, epsrel and restar.
c
c   on return
c
c     result real array of dimension numfun.
c            approximations to all components of the integral.
c     abserr real array of dimension numfun.
c            estimates of absolute errors.
c     neval  integer.
c            number of function evaluations used by dcuhre.
c     ifail  integer.
c            ifail = 0 for normal exit, when abserr(k) <=  epsabs or
c              abserr(k) <=  abs(result(k))*epsrel with maxpts or less
c              function evaluations for all values of k,
c              1 <= k <= numfun .
c            ifail = 1 if maxpts was too small for dcuhre
c              to obtain the required accuracy. in this case dcuhre
c              returns values of result with estimated absolute
c              errors abserr.
c            ifail = 2 if key is less than 0 or key greater than 4.
c            ifail = 3 if ndim is less than 2 or ndim greater than 15.
c            ifail = 4 if key = 1 and ndim not equal to 2.
c            ifail = 5 if key = 2 and ndim not equal to 3.
c            ifail = 6 if numfun is less than 1.
c            ifail = 7 if volume of region of integration is zero.
c            ifail = 8 if maxpts is less than 3*num.
c            ifail = 9 if maxpts is less than minpts.
c            ifail = 10 if epsabs < 0 and epsrel < 0.
c            ifail = 11 if nw is too small.
c            ifail = 12 if unlegal restar.
c     work   real array of dimension nw.
c            used as working storage.
c            work(nw) = nsub, the number of subregions in the data
c            structure.
c            let wrksub=(nw-1-17*numfun*mdiv)/(2*ndim+2*numfun+2)
c            work(1),...,work(numfun*wrksub) contain
c              the estimated components of the integrals over the
c              subregions.
c            work(numfun*wrksub+1),...,work(2*numfun*wrksub) contain
c              the estimated errors over the subregions.
c            work(2*numfun*wrksub+1),...,work(2*numfun*wrksub+ndim*
c              wrksub) contain the centers of the subregions.
c            work(2*numfun*wrksub+ndim*wrksub+1),...,work((2*numfun+
c              ndim)*wrksub+ndim*wrksub) contain subregion half widths.
c            work(2*numfun*wrksub+2*ndim*wrksub+1),...,work(2*numfun*
c              wrksub+2*ndim*wrksub+wrksub) contain the greatest errors
c              in each subregion.
c            work((2*numfun+2*ndim+1)*wrksub+1),...,work((2*numfun+
c              2*ndim+1)*wrksub+wrksub) contain the direction of
c              subdivision in each subregion.
c            work(2*(ndim+numfun+1)*wrksub),...,work(2*(ndim+numfun+1)*
c              wrksub+ 17*mdiv*numfun) is used as temporary
c              storage in dadhre.
c
c
c        dcuhre example test program
c
c 
c   dtest1 is a simple test driver for dcuhre. 
c 
c   output produced on a sun 3/50. 
c 
c       dcuhre test results 
c 
c    ftest calls = 3549, ifail =  0 
c   n   estimated error    integral 
c   1       0.000000       0.138508 
c   2       0.000000       0.063695 
c   3       0.000009       0.058618 
c   4       0.000000       0.054070 
c   5       0.000000       0.050056 
c   6       0.000000       0.046546 
c
c     program dtest1
c     external ftest
c     integer key, n, nf, ndim, mincls, maxcls, ifail, neval, nw
c     parameter (ndim = 5, nw = 5000, nf = ndim+1)
c     double precision a(ndim), b(ndim), wrkstr(nw)
c     double precision absest(nf), finest(nf), absreq, relreq
c     do 10 n = 1,ndim
c        a(n) = 0
c        b(n) = 1
c  10 continue
c     mincls = 0
c     maxcls = 10000
c     key = 0
c     absreq = 0
c     relreq = 1e-3
c     call dcuhre(ndim, nf, a, b, mincls, maxcls, ftest, absreq, relreq,
c    * key, nw, 0, finest, absest, neval, ifail, wrkstr)
c     print 9999, neval, ifail
c9999 format (8x, 'dcuhre test results', //'     ftest calls = ', i4,
c    * ', ifail = ', i2, /'    n   estimated error    integral')
c     do 20 n = 1,nf
c        print 9998, n, absest(n), finest(n)
c9998    format (3x, i2, 2f15.6)
c  20 continue
c     end
c     subroutine ftest(ndim, z, nfun, f)
c     integer n, ndim, nfun
c     double precision z(ndim), f(nfun), sum
c     sum = 0
c     do 10 n = 1,ndim
c        sum = sum + n*z(n)**2
c  10 continue
c     f(1) = exp(-sum/2)
c     do 20 n = 1,ndim
c        f(n+1) = z(n)*f(1)
c  20 continue
c     end
c
c***long description
c
c   the information for each subregion is contained in the
c   data structure work.
c   when passed on to dadhre, work is split into eight
c   arrays values, errors, centrs, hwidts, greate, dir,
c   oldres and work.
c   values contains the estimated values of the integrals.
c   errors contains the estimated errors.
c   centrs contains the centers of the subregions.
c   hwidts contains the half widths of the subregions.
c   greate contains the greatest estimated error for each subregion.
c   dir    contains the directions for further subdivision.
c   oldres and work are used as work arrays in dadhre.
c
c   the data structures for the subregions are in dadhre organized
c   as a heap, and the size of greate(i) defines the position of
c   region i in the heap. the heap is maintained by the program
c   dtrhre.
c
c   the subroutine dadhre is written for efficient execution on shared
c   memory parallel computer. on a computer with nproc processors we will
c   in each subdivision step divide mdiv regions, where mdiv is
c   chosen such that mod(2*mdiv,nproc) = 0, in totally 2*mdiv new regions.
c   each processor will then compute estimates of the integrals and errors
c   over 2*mdiv/nproc subregions in each subdivision step.
c   the subroutine for estimating the integral and the error over
c   each subregion, drlhre, uses work2 as a work array.
c   we must make sure that each processor writes its results to
c   separate parts of the memory, and therefore the sizes of work and
c   work2 are functions of mdiv.
c   in order to achieve parallel processing of subregions, compiler
c   directives should be placed in front of the do 200
c   loop in dadhre on machines like alliant and cray.
c
c***references
c   j.berntsen, t.o.espelid and a.genz, an adaptive algorithm
c   for the approximate calculation of multiple integrals,
c   to be published.
c
c   j.berntsen, t.o.espelid and a.genz, dcuhre: an adaptive
c   multidimensional integration routine for a vector of
c   integrals, to be published.
c
c***routines called dchhre,dadhre
c***end prologue dcuhre
c
c   global variables.
c
      external funsub
      integer ndim,numfun,minpts,maxpts,key,nw,restar
      integer neval,ifail
      double precision a(ndim),b(ndim),epsabs,epsrel
      double precision result(numfun),abserr(numfun),work(nw)
c
c   local variables.
c
c   mdiv   integer.
c          mdiv is the number of subregions that are divided in
c          each subdivision step in dadhre.
c          mdiv is chosen default to 1.
c          for efficient execution on parallel computers
c          with nproc processors mdiv should be set equal to
c          the smallest integer such that mod(2*mdiv,nproc) = 0.
c   maxdim integer.
c          the maximum allowed value of ndim.
c   maxwt  integer. the maximum number of weights used by the
c          integration rule.
c   wtleng integer.
c          the number of generators used by the selected rule.
c   work2  real work space. the length
c          depends on the parameters mdiv,maxdim and maxwt.
c   maxsub integer.
c          the maximum allowed number of subdivisions
c          for the given values of key, ndim and maxpts.
c   minsub integer.
c          the minimum allowed number of subregions for the given
c          values of minpts, key and ndim.
c   wrksub integer.
c          the maximum allowed number of subregions as a function
c          of nw, numfun, ndim and mdiv. this determines the length
c          of the main work arrays.
c   num    integer. the number of integrand evaluations needed
c          over each subregion.
c
      integer mdiv,maxwt,wtleng,maxdim,lenw2,maxsub,minsub
      integer num,nsub,lenw,keyf
      parameter (mdiv=1)
      parameter (maxdim=15)
      parameter (maxwt=14)
      parameter (lenw2=2*mdiv*maxdim* (maxwt+1)+12*maxwt+2*maxdim)
      integer wrksub,i1,i2,i3,i4,i5,i6,i7,i8,k1,k2,k3,k4,k5,k6,k7,k8
      double precision work2(lenw2)
c
c***first executable statement dcuhre
c
c   compute num, wtleng, maxsub and minsub,
c   and check the input parameters.
c
      call dchhre(maxdim,ndim,numfun,mdiv,a,b,minpts,maxpts,epsabs,
     +            epsrel,key,nw,restar,num,maxsub,minsub,keyf,
     +            ifail,wtleng)
      wrksub = (nw - 1 - 17*mdiv*numfun)/(2*ndim + 2*numfun + 2)
      if (ifail.ne.0) then
          go to 999
      end if
c
c   split up the work space.
c
      i1 = 1
      i2 = i1 + wrksub*numfun
      i3 = i2 + wrksub*numfun
      i4 = i3 + wrksub*ndim
      i5 = i4 + wrksub*ndim
      i6 = i5 + wrksub
      i7 = i6 + wrksub
      i8 = i7 + numfun*mdiv
      k1 = 1
      k2 = k1 + 2*mdiv*wtleng*ndim
      k3 = k2 + wtleng*5
      k4 = k3 + wtleng
      k5 = k4 + ndim
      k6 = k5 + ndim
      k7 = k6 + 2*mdiv*ndim
      k8 = k7 + 3*wtleng
c
c   on restart runs the number of subregions from the
c   previous call is assigned to nsub.
c
      if (restar.eq.1) then
          nsub = work(nw)
      end if
c
c   compute the size of the temporary work space needed in dadhre.
c
      lenw = 16*mdiv*numfun
      call dadhre(ndim,numfun,mdiv,a,b,minsub,maxsub,funsub,epsabs,
     +            epsrel,keyf,restar,num,lenw,wtleng,
     +            result,abserr,neval,nsub,ifail,work(i1),work(i2),
     +            work(i3),work(i4),work(i5),work(i6),work(i7),work(i8),
     +            work2(k1),work2(k2),work2(k3),work2(k4),work2(k5),
     +            work2(k6),work2(k7),work2(k8))
      work(nw) = nsub
999   return
c
c***end dcuhre
c
      end
      subroutine dchhre(maxdim,ndim,numfun,mdiv,a,b,minpts,maxpts,
     +                  epsabs,epsrel,key,nw,restar,num,maxsub,minsub,
     +                  keyf,ifail,wtleng)
c***begin prologue dchhre
c***purpose  dchhre checks the validity of the
c            input parameters to dcuhre.
c***description
c            dchhre computes num, maxsub, minsub, keyf, wtleng and 
c            ifail as functions of the input parameters to dcuhre,
c            and checks the validity of the input parameters to dcuhre.
c
c   on entry
c
c     maxdim integer.
c            the maximum allowed number of dimensions.
c     ndim   integer.
c            number of variables. 1 < ndim <= maxdim.
c     numfun integer.
c            number of components of the integral.
c     mdiv   integer.
c            mdiv is the number of subregions that are divided in
c            each subdivision step in dadhre.
c            mdiv is chosen default to 1.
c            for efficient execution on parallel computers
c            with nproc processors mdiv should be set equal to
c            the smallest integer such that mod(2*mdiv,nproc) = 0.
c     a      real array of dimension ndim.
c            lower limits of integration.
c     b      real array of dimension ndim.
c            upper limits of integration.
c     minpts integer.
c            minimum number of function evaluations.
c     maxpts integer.
c            maximum number of function evaluations.
c            the number of function evaluations over each subregion
c            is num.
c            if (key = 0 or key = 1) and (ndim = 2) then
c              num = 65
c            elseif (key = 0 or key = 2) and (ndim = 3) then
c              num = 127
c            elseif (key = 0 and ndim > 3) or (key = 3) then
c              num = 1 + 4*2*ndim + 2*ndim*(ndim-1) + 4*ndim*(ndim-1) +
c                    4*ndim*(ndim-1)*(ndim-2)/3 + 2**ndim
c            elseif (key = 4) then
c              num = 1 + 3*2*ndim + 2*ndim*(ndim-1) + 2**ndim
c            maxpts >= 3*num and maxpts >= minpts
c     epsabs real.
c            requested absolute error.
c     epsrel real.
c            requested relative error.
c     key    integer.
c            key to selected local integration rule.
c            key = 0 is the default value.
c                  for ndim = 2 the degree 13 rule is selected.
c                  for ndim = 3 the degree 11 rule is selected.
c                  for ndim > 3 the degree  9 rule is selected.
c            key = 1 gives the user the 2 dimensional degree 13
c                  integration rule that uses 65 evaluation points.
c            key = 2 gives the user the 3 dimensional degree 11
c                  integration rule that uses 127 evaluation points.
c            key = 3 gives the user the degree 9 integration rule.
c            key = 4 gives the user the degree 7 integration rule.
c                  this is the recommended rule for problems that
c                  require great adaptivity.
c     nw     integer.
c            defines the length of the working array work.
c            let maxsub denote the maximum allowed number of subregions
c            for the given values of maxpts, key and ndim.
c            maxsub = (maxpts-num)/(2*num) + 1
c            nw should be greater or equal to
c            maxsub*(2*ndim+2*numfun+2) + 17*numfun + 1
c            for efficient execution on parallel computers
c            nw should be chosen greater or equal to
c            maxsub*(2*ndim+2*numfun+2) + 17*numfun*mdiv + 1
c            where mdiv is the number of subregions that are divided in
c            each subdivision step.
c            mdiv is default set internally in dcuhre equal to 1.
c            for efficient execution on parallel computers
c            with nproc processors mdiv should be set equal to
c            the smallest integer such that mod(2*mdiv,nproc) = 0.
c     restar integer.
c            if restar = 0, this is the first attempt to compute
c            the integral.
c            if restar = 1, then we restart a previous attempt.
c
c   on return
c
c     num    integer.
c            the number of function evaluations over each subregion.
c     maxsub integer.
c            the maximum allowed number of subregions for the
c            given values of maxpts, key and ndim.
c     minsub integer.
c            the minimum allowed number of subregions for the given
c            values of minpts, key and ndim.
c     ifail  integer.
c            ifail = 0 for normal exit.
c            ifail = 2 if key is less than 0 or key greater than 4.
c            ifail = 3 if ndim is less than 2 or ndim greater than
c                      maxdim.
c            ifail = 4 if key = 1 and ndim not equal to 2.
c            ifail = 5 if key = 2 and ndim not equal to 3.
c            ifail = 6 if numfun less than 1.
c            ifail = 7 if volume of region of integration is zero.
c            ifail = 8 if maxpts is less than 3*num.
c            ifail = 9 if maxpts is less than minpts.
c            ifail = 10 if epsabs < 0 and epsrel < 0.
c            ifail = 11 if nw is too small.
c            ifail = 12 if unlegal restar.
c     keyf   integer.
c            key to selected integration rule.
c     wtleng integer.
c            the number of generators of the chosen integration rule.
c
c***routines called-none
c***end prologue dchhre
c
c   global variables.
c
      integer ndim,numfun,mdiv,minpts,maxpts,key,nw,minsub,maxsub
      integer restar,num,keyf,ifail,maxdim,wtleng
      double precision a(ndim),b(ndim),epsabs,epsrel
c
c   local variables.
c
      integer limit,j
c
c***first executable statement dchhre
c
      ifail = 0
c
c   check on legal key.
c
      if (key.lt.0 .or. key.gt.4) then
          ifail = 2
          go to 999
      end if
c
c   check on legal ndim.
c
      if (ndim.lt.2 .or. ndim.gt.maxdim) then
          ifail = 3
          go to 999
      end if
c
c   for key = 1, ndim must be equal to 2.
c
      if (key.eq.1 .and. ndim.ne.2) then
          ifail = 4
          go to 999
      end if
c
c   for key = 2, ndim must be equal to 3.
c
      if (key.eq.2 .and. ndim.ne.3) then
          ifail = 5
          go to 999
      end if
c
c   for key = 0, we point at the selected integration rule.
c
      if (key.eq.0) then
          if (ndim.eq.2) then
              keyf = 1
          else if (ndim.eq.3) then
              keyf = 2
          else
              keyf = 3
          endif
      else
          keyf = key
      endif
c
c   compute num and wtleng as a function of keyf and ndim.
c
      if (keyf.eq.1) then
          num = 65
          wtleng = 14
      else if (keyf.eq.2) then
          num = 127
          wtleng = 13
      else if (keyf.eq.3) then
          num = 1 + 4*2*ndim + 2*ndim* (ndim-1) + 4*ndim* (ndim-1) +
     +          4*ndim* (ndim-1)* (ndim-2)/3 + 2**ndim
          wtleng = 9
          if (ndim.eq.2) wtleng = 8
      else if (keyf.eq.4) then
          num = 1 + 3*2*ndim + 2*ndim* (ndim-1) + 2**ndim
          wtleng = 6
      end if
c
c   compute maxsub.
c
      maxsub = (maxpts-num)/ (2*num) + 1
c
c   compute minsub.
c
      minsub = (minpts-num)/ (2*num) + 1
      if (mod(minpts-num,2*num).ne.0) then
          minsub = minsub + 1
      end if
      minsub = max(2,minsub)
c
c   check on positive numfun.
c
      if (numfun.lt.1) then
          ifail = 6
          go to 999
      end if
c
c   check on legal upper and lower limits of integration.
c
      do 10 j = 1,ndim
          if (a(j)-b(j).eq.0) then
              ifail = 7
              go to 999
          end if
10    continue
c
c   check on maxpts < 3*num.
c
      if (maxpts.lt.3*num) then
          ifail = 8
          go to 999
      end if
c
c   check on maxpts >= minpts.
c
      if (maxpts.lt.minpts) then
          ifail = 9
          go to 999
      end if
c
c   check on legal accuracy requests.
c
      if (epsabs.lt.0 .and. epsrel.lt.0) then
          ifail = 10
          go to 999
      end if
c
c   check on big enough double precision workspace.
c
      limit = maxsub* (2*ndim+2*numfun+2) + 17*mdiv*numfun + 1
      if (nw.lt.limit) then
          ifail = 11
          go to 999
      end if
c
c    check on legal restar.
c
      if (restar.ne.0 .and. restar.ne.1) then
          ifail = 12
          go to 999
      end if
999   return
c
c***end dchhre
c
      end
      subroutine dadhre(ndim,numfun,mdiv,a,b,minsub,maxsub,funsub,
     +                  epsabs,epsrel,key,restar,num,lenw,wtleng,
     +                  result,abserr,neval,nsub,ifail,values,
     +                  errors,centrs,hwidts,greate,dir,oldres,work,g,w,
     +                  rulpts,center,hwidth,x,scales,norms)
c***begin prologue dadhre
c***keywords automatic multidimensional integrator,
c            n-dimensional hyper-rectangles,
c            general purpose, global adaptive
c***purpose  the routine calculates an approximation to a given
c            vector of definite integrals, i, over a hyper-rectangular
c            region hopefully satisfying for each component of i the
c            following claim for accuracy:
c            abs(i(k)-result(k)).le.max(epsabs,epsrel*abs(i(k)))
c***description computation of integrals over hyper-rectangular
c            regions.
c            dadhre repeatedly subdivides the region
c            of integration and estimates the integrals and the
c            errors over the subregions with  greatest
c            estimated errors until the error request
c            is met or maxsub subregions are stored.
c            the regions are divided in two equally sized parts along
c            the direction with greatest absolute fourth divided
c            difference.
c
c   on entry
c
c     ndim   integer.
c            number of variables. 1 < ndim <= maxdim.
c     numfun integer.
c            number of components of the integral.
c     mdiv   integer.
c            defines the number of new subregions that are divided
c            in each subdivision step.
c     a      real array of dimension ndim.
c            lower limits of integration.
c     b      real array of dimension ndim.
c            upper limits of integration.
c     minsub integer.
c            the computations proceed until there are at least
c            minsub subregions in the data structure.
c     maxsub integer.
c            the computations proceed until there are at most
c            maxsub subregions in the data structure.
c
c     funsub externally declared subroutine for computing
c            all components of the integrand in the given
c            evaluation point.
c            it must have parameters (ndim,x,numfun,funvls)
c            input parameters:
c              ndim   integer that defines the dimension of the
c                     integral.
c              x      real array of dimension ndim
c                     that defines the evaluation point.
c              numfun integer that defines the number of
c                     components of i.
c            output parameter:
c              funvls real array of dimension numfun
c                     that defines numfun components of the integrand.
c
c     epsabs real.
c            requested absolute error.
c     epsrel real.
c            requested relative error.
c     key    integer.
c            key to selected local integration rule.
c            key = 0 is the default value.
c                  for ndim = 2 the degree 13 rule is selected.
c                  for ndim = 3 the degree 11 rule is selected.
c                  for ndim > 3 the degree  9 rule is selected.
c            key = 1 gives the user the 2 dimensional degree 13
c                  integration rule that uses 65 evaluation points.
c            key = 2 gives the user the 3 dimensional degree 11
c                  integration rule that uses 127 evaluation points.
c            key = 3 gives the user the degree 9 integration rule.
c            key = 4 gives the user the degree 7 integration rule.
c                  this is the recommended rule for problems that
c                  require great adaptivity.
c     restar integer.
c            if restar = 0, this is the first attempt to compute
c            the integral.
c            if restar = 1, then we restart a previous attempt.
c            (in this case the output parameters from dadhre
c            must not be changed since the last
c            exit from dadhre.)
c     num    integer.
c            the number of function evaluations over each subregion.
c     lenw   integer.
c            defines the length of the working array work.
c            lenw should be greater or equal to
c            16*mdiv*numfun.
c     wtleng integer.
c            the number of weights in the basic integration rule.
c     nsub   integer.
c            if restar = 1, then nsub must specify the number
c            of subregions stored in the previous call to dadhre.
c
c   on return
c
c     result real array of dimension numfun.
c            approximations to all components of the integral.
c     abserr real array of dimension numfun.
c            estimates of absolute accuracies.
c     neval  integer.
c            number of function evaluations used by dadhre.
c     nsub   integer.
c            number of stored subregions.
c     ifail  integer.
c            ifail = 0 for normal exit, when abserr(k) <=  epsabs or
c              abserr(k) <=  abs(result(k))*epsrel with maxsub or less
c              subregions processed for all values of k,
c              1 <=  k <=  numfun.
c            ifail = 1 if maxsub was too small for dadhre
c              to obtain the required accuracy. in this case dadhre
c              returns values of result with estimated absolute
c              accuracies abserr.
c     values real array of dimension (numfun,maxsub).
c            used to store estimated values of the integrals
c            over the subregions.
c     errors real array of dimension (numfun,maxsub).
c            used to store the corresponding estimated errors.
c     centrs real array of dimension (ndim,maxsub).
c            used to store the centers of the stored subregions.
c     hwidts real array of dimension (ndim,maxsub).
c            used to store the half widths of the stored subregions.
c     greate real array of dimension maxsub.
c            used to store the greatest estimated errors in
c            all subregions.
c     dir    real array of dimension maxsub.
c            dir is used to store the directions for
c            further subdivision.
c     oldres real array of dimension (numfun,mdiv).
c            used to store old estimates of the integrals over the
c            subregions.
c     work   real array of dimension lenw.
c            used  in drlhre and dtrhre.
c     g      real array of dimension (ndim,wtleng,2*mdiv).
c            the fully symmetric sum generators for the rules.
c            g(1,j,1),...,g(ndim,j,1) are the generators for the
c            points associated with the jth weights.
c            when mdiv subregions are divided in 2*mdiv
c            subregions, the subregions may be processed on different
c            processors and we must make a copy of the generators
c            for each processor.
c     w      real array of dimension (5,wtleng).
c            the weights for the basic and null rules.
c            w(1,1), ..., w(1,wtleng) are weights for the basic rule.
c            w(i,1), ..., w(i,wtleng) , for i > 1 are null rule weights.
c     rulpts real array of dimension wtleng.
c            work array used in dinhre.
c     center real array of dimension ndim.
c            work array used in dtrhre.
c     hwidth real array of dimension ndim.
c            work array used in dtrhre.
c     x      real array of dimension (ndim,2*mdiv).
c            work array used in drlhre.
c     scales real array of dimension (3,wtleng).
c            work array used by dinhre and drlhre.
c     norms  real array of dimension (3,wtleng).
c            work array used by dinhre and drlhre.
c
c***references
c
c   p. van dooren and l. de ridder, algorithm 6, an adaptive algorithm
c   for numerical integration over an n-dimensional cube, j.comput.appl.
c   math. 2(1976)207-217.
c
c   a.c.genz and a.a.malik, algorithm 019. remarks on algorithm 006:
c   an adaptive algorithm for numerical integration over an
c   n-dimensional rectangular region,j.comput.appl.math. 6(1980)295-302.
c
c***routines called dtrhre,dinhre,drlhre
c***end prologue dadhre
c
c   global variables.
c
      external funsub
      integer ndim,numfun,mdiv,minsub,maxsub,key,lenw,restar
      integer num,neval,nsub,ifail,wtleng
      double precision a(ndim),b(ndim),epsabs,epsrel
      double precision result(numfun),abserr(numfun)
      double precision values(numfun,maxsub),errors(numfun,maxsub)
      double precision centrs(ndim,maxsub)
      double precision hwidts(ndim,maxsub)
      double precision greate(maxsub),dir(maxsub)
      double precision oldres(numfun,mdiv)
      double precision work(lenw),rulpts(wtleng)
      double precision g(ndim,wtleng,2*mdiv),w(5,wtleng)
      double precision center(ndim),hwidth(ndim),x(ndim,2*mdiv)
      double precision scales(3,wtleng),norms(3,wtleng)
c
c   local variables.
c
c   intsgn is used to get correct sign on the integral.
c   sbrgns is the number of stored subregions.
c   ndiv   the number of subregions to be divided in each main step.
c   pointr pointer to the position in the data structure where
c          the new subregions are to be stored.
c   direct direction of subdivision.
c   errcof heuristic error coeff. defined in dinhre and used by drlhre
c          and dadhre.
c
      integer i,j,k
      integer intsgn,sbrgns
      integer l1
      integer ndiv,pointr,direct,index
      double precision oldcen,est1,est2,errcof(6)
c
c***first executable statement dadhre
c
c   get the correct sign on the integral.
c
      intsgn = 1
      do 10 j = 1,ndim
          if (b(j).lt.a(j)) then
              intsgn = - intsgn
          end if
10    continue
c
c   call dinhre to compute the weights and abscissas of
c   the function evaluation points.
c
      call dinhre(ndim,key,wtleng,w,g,errcof,rulpts,scales,norms)
c
c   if restar = 1, then this is a restart run.
c
      if (restar.eq.1) then
          sbrgns = nsub
          go to 110
      end if
c
c   initialize the sbrgns, centrs and hwidts.
c
      sbrgns = 1
      do 15 j = 1,ndim
          centrs(j,1) = (a(j)+b(j))/2
          hwidts(j,1) = abs(b(j)-a(j))/2
15    continue
c
c   initialize result, abserr and neval.
c
      do 20 j = 1,numfun
          result(j) = 0
          abserr(j) = 0
20    continue
      neval = 0
c
c   apply drlhre over the whole region.
c
      call drlhre(ndim,centrs(1,1),hwidts(1,1),wtleng,g,w,errcof,numfun,
     +            funsub,scales,norms,x,work,values(1,1),errors(1,1),
     +            dir(1))
      neval = neval + num
c
c   add the computed values to result and abserr.
c
      do 55 j = 1,numfun
          result(j) = result(j) + values(j,1)
55    continue
      do 65 j = 1,numfun
          abserr(j) = abserr(j) + errors(j,1)
65    continue
c
c   store results in heap.
c
      index = 1
      call dtrhre(2,ndim,numfun,index,values,errors,centrs,hwidts,
     +            greate,work(1),work(numfun+1),center,hwidth,dir)
c
c***end initialisation.
c
c***begin loop while the error is too great,
c   and sbrgns+1 is less than maxsub.
c
110   if (sbrgns+1.le.maxsub) then
c
c   if we are allowed to divide further,
c   prepare to apply basic rule over each half of the
c   ndiv subregions with greatest errors.
c   if maxsub is great enough, ndiv = mdiv
c
          if (mdiv.gt.1) then
              ndiv = maxsub - sbrgns
              ndiv = min(ndiv,mdiv,sbrgns)
          else
              ndiv = 1
          end if
c
c   divide the ndiv subregions in two halves, and compute
c   integral and error over each half.
c
          do 150 i = 1,ndiv
              pointr = sbrgns + ndiv + 1 - i
c
c   adjust result and abserr.
c
              do 115 j = 1,numfun
                  result(j) = result(j) - values(j,1)
                  abserr(j) = abserr(j) - errors(j,1)
115           continue
c
c   compute first half region.
c
              do 120 j = 1,ndim
                  centrs(j,pointr) = centrs(j,1)
                  hwidts(j,pointr) = hwidts(j,1)
120           continue
              direct = dir(1)
              dir(pointr) = direct
              hwidts(direct,pointr) = hwidts(direct,1)/2
              oldcen = centrs(direct,1)
              centrs(direct,pointr) = oldcen - hwidts(direct,pointr)
c
c   save the computed values of the integrals.
c
              do 125 j = 1,numfun
                  oldres(j,ndiv-i+1) = values(j,1)
125           continue
c
c   adjust the heap.
c
              call dtrhre(1,ndim,numfun,sbrgns,values,errors,centrs,
     +                    hwidts,greate,work(1),work(numfun+1),center,
     +                    hwidth,dir)
c
c   compute second half region.
c
              do 130 j = 1,ndim
                  centrs(j,pointr-1) = centrs(j,pointr)
                  hwidts(j,pointr-1) = hwidts(j,pointr)
130           continue
              centrs(direct,pointr-1) = oldcen + hwidts(direct,pointr)
              hwidts(direct,pointr-1) = hwidts(direct,pointr)
              dir(pointr-1) = direct
150       continue
c
c   make copies of the generators for each processor.
c
          do 190 i = 2,2*ndiv
              do 190 j = 1,ndim
                  do 190 k = 1,wtleng
                      g(j,k,i) = g(j,k,1)
190       continue
c
c   apply basic rule.
c
cvd$l cncall
          do 200 i = 1,2*ndiv
              index = sbrgns + i
              l1 = 1 + (i-1)*8*numfun
              call drlhre(ndim,centrs(1,index),hwidts(1,index),wtleng,
     +                    g(1,1,i),w,errcof,numfun,funsub,scales,norms,
     +                    x(1,i),work(l1),values(1,index),
     +                    errors(1,index),dir(index))
200       continue
          neval = neval + 2*ndiv*num
c
c   add new contributions to result.
c
          do 220 i = 1,2*ndiv
              do 210 j = 1,numfun
                  result(j) = result(j) + values(j,sbrgns+i)
210           continue
220       continue
c
c   check consistency of results and if necessary adjust
c   the estimated errors.
c
          do 240 i = 1,ndiv
              greate(sbrgns+2*i-1) = 0
              greate(sbrgns+2*i) = 0
              do 230 j = 1,numfun
                  est1 = abs(oldres(j,i)- (values(j,
     +                   sbrgns+2*i-1)+values(j,sbrgns+2*i)))
                  est2 = errors(j,sbrgns+2*i-1) + errors(j,sbrgns+2*i)
                  if (est2.gt.0) then
                      errors(j,sbrgns+2*i-1) = errors(j,sbrgns+2*i-1)*
     +                  (1+errcof(5)*est1/est2)
                      errors(j,sbrgns+2*i) = errors(j,sbrgns+2*i)*
     +                                       (1+errcof(5)*est1/est2)
                  end if
                  errors(j,sbrgns+2*i-1) = errors(j,sbrgns+2*i-1) +
     +                                     errcof(6)*est1
                  errors(j,sbrgns+2*i) = errors(j,sbrgns+2*i) +
     +                                   errcof(6)*est1
                  if (errors(j,sbrgns+2*i-1).gt.
     +                greate(sbrgns+2*i-1)) then
                      greate(sbrgns+2*i-1) = errors(j,sbrgns+2*i-1)
                  end if
                  if (errors(j,sbrgns+2*i).gt.greate(sbrgns+2*i)) then
                      greate(sbrgns+2*i) = errors(j,sbrgns+2*i)
                  end if
                  abserr(j) = abserr(j) + errors(j,sbrgns+2*i-1) +
     +                        errors(j,sbrgns+2*i)
230           continue
240       continue
c
c   store results in heap.
c
          do 250 i = 1,2*ndiv
              index = sbrgns + i
              call dtrhre(2,ndim,numfun,index,values,errors,centrs,
     +                    hwidts,greate,work(1),work(numfun+1),center,
     +                    hwidth,dir)
250       continue
          sbrgns = sbrgns + 2*ndiv
c
c   check for termination.
c
          if (sbrgns.lt.minsub) then
              go to 110
          end if
          do 255 j = 1,numfun
              if (abserr(j).gt.epsrel*abs(result(j)) .and.
     +            abserr(j).gt.epsabs) then
                  go to 110
              end if
255       continue
          ifail = 0
          go to 499
c
c   else we did not succeed with the
c   given value of maxsub.
c
      else
          ifail = 1
      end if
c
c   compute more accurate values of result and abserr.
c
499   continue
      do 500 j = 1,numfun
          result(j) = 0
          abserr(j) = 0
500   continue
      do 510 i = 1,sbrgns
          do 505 j = 1,numfun
              result(j) = result(j) + values(j,i)
              abserr(j) = abserr(j) + errors(j,i)
505       continue
510   continue
c
c   compute correct sign on the integral.
c
      do 600 j = 1,numfun
          result(j) = result(j)*intsgn
600   continue
      nsub = sbrgns
      return
c
c***end dadhre
c
      end
      subroutine dinhre(ndim,key,wtleng,w,g,errcof,rulpts,scales,norms)
c***begin prologue dinhre
c***purpose dinhre computes abscissas and weights of the integration
c            rule and the null rules to be used in error estimation.
c            these are computed as functions of ndim and key.
c***description dinhre will for given values of ndim and key compute or
c            select the correct values of the abscissas and
c            corresponding weights for different
c            integration rules and null rules and assign them to
c            g and w.
c            the heuristic error coefficients errcof
c            will be computed as a function of key.
c            scaling factors scales and normalization factors norms
c            used in the error estimation are computed.
c
c
c   on entry
c
c     ndim   integer.
c            number of variables.
c     key    integer.
c            key to selected local integration rule.
c     wtleng integer.
c            the number of weights in each of the rules.
c
c   on return
c
c     w      real array of dimension (5,wtleng).
c            the weights for the basic and null rules.
c            w(1,1), ...,w(1,wtleng) are weights for the basic rule.
c            w(i,1), ...,w(i,wtleng), for i > 1 are null rule weights.
c     g      real array of dimension (ndim,wtleng).
c            the fully symmetric sum generators for the rules.
c            g(1,j),...,g(ndim,j) are the generators for the points
c            associated with the the jth weights.
c     errcof real array of dimension 6.
c            heuristic error coefficients that are used in the
c            error estimation in basrul.
c            it is assumed that the error is computed using:
c             if (n1*errcof(1) < n2 and n2*errcof(2) < n3)
c               then error = errcof(3)*n1
c               else error = errcof(4)*max(n1, n2, n3)
c             error = error + ep*(errcof(5)*error/(es+error)+errcof(6))
c            where n1-n3 are the null rules, ep is the error for
c            the parent
c            subregion and es is the error for the sibling subregion.
c     rulpts real array of dimension wtleng.
c            a work array containing the number of points produced by
c            each generator of the selected rule.
c     scales real array of dimension (3,wtleng).
c            scaling factors used to construct new null rules,
c            n1, n2 and n3,
c            based on a linear combination of two successive null rules
c            in the sequence of null rules.
c     norms  real array of dimension (3,wtleng).
c            2**ndim/(1-norm of the null rule constructed by each of
c            the scaling factors.)
c
c***routines called  d132re,d113re,d07hre,d09hre
c***end prologue dinhre
c
c   global variables.
c
      integer ndim,key,wtleng
      double precision g(ndim,wtleng),w(5,wtleng),errcof(6)
      double precision rulpts(wtleng),scales(3,wtleng)
      double precision norms(3,wtleng)
c
c   local variables.
c
      integer i,j,k
      double precision we(14)
c
c***first executable statement dinhre
c
c   compute w, g and errcof.
c
      if (key.eq.1) then
          call d132re(wtleng,w,g,errcof,rulpts)
      else if (key.eq.2) then
          call d113re(wtleng,w,g,errcof,rulpts)
      else if (key.eq.3) then
          call d09hre(ndim,wtleng,w,g,errcof,rulpts)
      else if (key.eq.4) then
          call d07hre(ndim,wtleng,w,g,errcof,rulpts)
      end if
c
c   compute scales and norms.
c
      do 100 k = 1,3
          do 50 i = 1,wtleng
              if (w(k+1,i).ne.0) then
                  scales(k,i) = - w(k+2,i)/w(k+1,i)
              else
                  scales(k,i) = 100
              end if
              do 30 j = 1,wtleng
                  we(j) = w(k+2,j) + scales(k,i)*w(k+1,j)
30            continue
              norms(k,i) = 0
              do 40 j = 1,wtleng
                  norms(k,i) = norms(k,i) + rulpts(j)*abs(we(j))
40            continue
              norms(k,i) = 2**ndim/norms(k,i)
50        continue
100   continue
      return
c
c***end dinhre
c
      end
      subroutine d132re(wtleng,w,g,errcof,rulpts)
c***begin prologue d132re
c***author   jarle berntsen, edb-senteret,
c            university of bergen, thormohlens gt. 55,
c            n-5008 bergen, norway
c***purpose d132re computes abscissas and weights of a 2 dimensional
c            integration rule of degree 13.
c            two null rules of degree 11, one null rule of degree 9
c            and one null rule of degree 7 to be used in error
c            estimation are also computed.
c ***description d132re will select the correct values of the abscissas
c            and corresponding weights for different
c            integration rules and null rules and assign them to
c            g and w. the heuristic error coefficients errcof
c            will also be assigned.
c
c
c   on entry
c
c     wtleng integer.
c            the number of weights in each of the rules.
c
c   on return
c
c     w      real array of dimension (5,wtleng).
c            the weights for the basic and null rules.
c            w(1,1),...,w(1,wtleng) are weights for the basic rule.
c            w(i,1),...,w(i,wtleng), for i > 1 are null rule weights.
c     g      real array of dimension (ndim,wtleng).
c            the fully symmetric sum generators for the rules.
c            g(1,j),...,g(ndim,j) are the generators for the points
c            associated with the the jth weights.
c     errcof real array of dimension 6.
c            heuristic error coefficients that are used in the
c            error estimation in basrul.
c     rulpts real array of dimension wtleng.
c            the number of points produced by each generator.
c***references s.eriksen,
c              thesis of the degree cand.scient, dept. of informatics,
c              univ. of bergen,norway, 1984.
c
c***routines called-none
c***end prologue d132re
c
c   global variables
c
      integer wtleng
      double precision w(5,wtleng),g(2,wtleng),errcof(6)
      double precision rulpts(wtleng)
c
c   local variables.
c
      integer i,j
      double precision dim2g(16)
      double precision dim2w(14,5)
c
      data (dim2g(i),i=1,16)/0.2517129343453109d+00,
     +     0.7013933644534266d+00,0.9590960631619962d+00,
     +     0.9956010478552127d+00,0.5000000000000000d+00,
     +     0.1594544658297559d+00,0.3808991135940188d+00,
     +     0.6582769255267192d+00,0.8761473165029315d+00,
     +     0.9982431840531980d+00,0.9790222658168462d+00,
     +     0.6492284325645389d+00,0.8727421201131239d+00,
     +     0.3582614645881228d+00,0.5666666666666666d+00,
     +     0.2077777777777778d+00/
c
      data (dim2w(i,1),i=1,14)/0.3379692360134460d-01,
     +     0.9508589607597761d-01,0.1176006468056962d+00,
     +     0.2657774586326950d-01,0.1701441770200640d-01,
     +     0.0000000000000000d+00,0.1626593098637410d-01,
     +     0.1344892658526199d+00,0.1328032165460149d+00,
     +     0.5637474769991870d-01,0.3908279081310500d-02,
     +     0.3012798777432150d-01,0.1030873234689166d+00,
     +     0.6250000000000000d-01/
c
      data (dim2w(i,2),i=1,14)/0.3213775489050763d+00,
     +     - .1767341636743844d+00,0.7347600537466072d-01,
     +     - .3638022004364754d-01,0.2125297922098712d-01,
     +     0.1460984204026913d+00,0.1747613286152099d-01,
     +     0.1444954045641582d+00,0.1307687976001325d-03,
     +     0.5380992313941161d-03,0.1042259576889814d-03,
     +     - .1401152865045733d-02,0.8041788181514763d-02,
     +     - .1420416552759383d+00/
c
      data (dim2w(i,3),i=1,14)/0.3372900883288987d+00,
     +     - .1644903060344491d+00,0.7707849911634622d-01,
     +     - .3804478358506310d-01,0.2223559940380806d-01,
     +     0.1480693879765931d+00,0.4467143702185814d-05,
     +     0.1508944767074130d+00,0.3647200107516215d-04,
     +     0.5777198999013880d-03,0.1041757313688177d-03,
     +     - .1452822267047819d-02,0.8338339968783705d-02,
     +     - .1472796329231960d+00/
c
      data (dim2w(i,4),i=1,14)/ - .8264123822525677d+00,
     +     0.3065838614094360d+00,0.2389292538329435d-02,
     +     - .1343024157997222d+00,0.8833366840533900d-01,
     +     0.0000000000000000d+00,0.9786283074168292d-03,
     +     - .1319227889147519d+00,0.7990012200150630d-02,
     +     0.3391747079760626d-02,0.2294915718283264d-02,
     +     - .1358584986119197d-01,0.4025866859057809d-01,
     +     0.3760268580063992d-02/
c
      data (dim2w(i,5),i=1,14)/0.6539094339575232d+00,
     +     - .2041614154424632d+00, - .1746981515794990d+00,
     +     0.3937939671417803d-01,0.6974520545933992d-02,
     +     0.0000000000000000d+00,0.6667702171778258d-02,
     +     0.5512960621544304d-01,0.5443846381278607d-01,
     +     0.2310903863953934d-01,0.1506937747477189d-01,
     +     - .6057021648901890d-01,0.4225737654686337d-01,
     +     0.2561989142123099d-01/
c
c***first executable statement d132re
c
c   assign values to w.
c
      do 10 i = 1,14
          do 10 j = 1,5
              w(j,i) = dim2w(i,j)
10    continue
c
c   assign values to g.
c
      do 20 i = 1,2
          do 20 j = 1,14
              g(i,j) = 0
20    continue
      g(1,2) = dim2g(1)
      g(1,3) = dim2g(2)
      g(1,4) = dim2g(3)
      g(1,5) = dim2g(4)
      g(1,6) = dim2g(5)
      g(1,7) = dim2g(6)
      g(2,7) = g(1,7)
      g(1,8) = dim2g(7)
      g(2,8) = g(1,8)
      g(1,9) = dim2g(8)
      g(2,9) = g(1,9)
      g(1,10) = dim2g(9)
      g(2,10) = g(1,10)
      g(1,11) = dim2g(10)
      g(2,11) = g(1,11)
      g(1,12) = dim2g(11)
      g(2,12) = dim2g(12)
      g(1,13) = dim2g(13)
      g(2,13) = dim2g(14)
      g(1,14) = dim2g(15)
      g(2,14) = dim2g(16)
c
c   assign values to rulpts.
c
      rulpts(1) = 1
      do 30 i = 2,11
          rulpts(i) = 4
30    continue
      rulpts(12) = 8
      rulpts(13) = 8
      rulpts(14) = 8
c
c   assign values to errcof.
c
      errcof(1) = 10
      errcof(2) = 10
      errcof(3) = 1.
      errcof(4) = 5.
      errcof(5) = 0.5
      errcof(6) = 0.25
c
c***end d132re
c
      return
      end
      subroutine d113re(wtleng,w,g,errcof,rulpts)
c***begin prologue d113re
c***author   jarle berntsen, edb-senteret,
c            university of bergen, thormohlens gt. 55,
c            n-5008 bergen, norway
c***purpose d113re computes abscissas and weights of a 3 dimensional
c            integration rule of degree 11.
c            two null rules of degree 9, one null rule of degree 7
c            and one null rule of degree 5 to be used in error
c            estimation are also computed.
c***description d113re will select the correct values of the abscissas
c            and corresponding weights for different
c            integration rules and null rules and assign them to g
c            and w.
c            the heuristic error coefficients errcof
c            will also be computed.
c
c
c   on entry
c
c     wtleng integer.
c            the number of weights in each of the rules.
c
c   on return
c
c     w      real array of dimension (5,wtleng).
c            the weights for the basic and null rules.
c            w(1,1),...,w(1,wtleng) are weights for the basic rule.
c            w(i,1),...,w(i,wtleng), for i > 1 are null rule weights.
c     g      real array of dimension (ndim,wtleng).
c            the fully symmetric sum generators for the rules.
c            g(1,j),...,g(ndim,j) are the generators for the points
c            associated with the the jth weights.
c     errcof real array of dimension 6.
c            heuristic error coefficients that are used in the
c            error estimation in basrul.
c     rulpts real array of dimension wtleng.
c            the number of points used by each generator.
c
c***references  j.berntsen, cautious adaptive numerical integration
c               over the 3-cube, reports in informatics 17, dept. of
c               inf.,univ. of bergen, norway, 1985.
c               j.berntsen and t.o.espelid, on the construction of
c               higher degree three-dimensional embedded integration
c               rules, siam j. numer. anal.,vol. 25,no. 1, pp.222-234,
c               1988.
c***routines called-none
c***end prologue d113re
c
c   global variables.
c
      integer wtleng
      double precision w(5,wtleng),g(3,wtleng),errcof(6)
      double precision rulpts(wtleng)
c
c   local variables.
c
      integer i,j
      double precision dim3g(14)
      double precision dim3w(13,5)
c
      data (dim3g(i),i=1,14)/0.1900000000000000d+00,
     +     0.5000000000000000d+00,0.7500000000000000d+00,
     +     0.8000000000000000d+00,0.9949999999999999d+00,
     +     0.9987344998351400d+00,0.7793703685672423d+00,
     +     0.9999698993088767d+00,0.7902637224771788d+00,
     +     0.4403396687650737d+00,0.4378478459006862d+00,
     +     0.9549373822794593d+00,0.9661093133630748d+00,
     +     0.4577105877763134d+00/
c
      data (dim3w(i,1),i=1,13)/0.7923078151105734d-02,
     +     0.6797177392788080d-01,0.1086986538805825d-02,
     +     0.1838633662212829d+00,0.3362119777829031d-01,
     +     0.1013751123334062d-01,0.1687648683985235d-02,
     +     0.1346468564512807d+00,0.1750145884600386d-02,
     +     0.7752336383837454d-01,0.2461864902770251d+00,
     +     0.6797944868483039d-01,0.1419962823300713d-01/
c
      data (dim3w(i,2),i=1,13)/0.1715006248224684d+01,
     +     - .3755893815889209d+00,0.1488632145140549d+00,
     +     - .2497046640620823d+00,0.1792501419135204d+00,
     +     0.3446126758973890d-02, - .5140483185555825d-02,
     +     0.6536017839876425d-02, - .6513454939229700d-03,
     +     - .6304672433547204d-02,0.1266959399788263d-01,
     +     - .5454241018647931d-02,0.4826995274768427d-02/
c
      data (dim3w(i,3),i=1,13)/0.1936014978949526d+01,
     +     - .3673449403754268d+00,0.2929778657898176d-01,
     +     - .1151883520260315d+00,0.5086658220872218d-01,
     +     0.4453911087786469d-01, - .2287828257125900d-01,
     +     0.2908926216345833d-01, - .2898884350669207d-02,
     +     - .2805963413307495d-01,0.5638741361145884d-01,
     +     - .2427469611942451d-01,0.2148307034182882d-01/
c
      data (dim3w(i,4),i=1,13)/0.5170828195605760d+00,
     +     0.1445269144914044d-01, - .3601489663995932d+00,
     +     0.3628307003418485d+00,0.7148802650872729d-02,
     +     - .9222852896022966d-01,0.1719339732471725d-01,
     +     - .1021416537460350d+00, - .7504397861080493d-02,
     +     0.1648362537726711d-01,0.5234610158469334d-01,
     +     0.1445432331613066d-01,0.3019236275367777d-02/
c
      data (dim3w(i,5),i=1,13)/0.2054404503818520d+01,
     +     0.1377759988490120d-01, - .5768062917904410d+00,
     +     0.3726835047700328d-01,0.6814878939777219d-02,
     +     0.5723169733851849d-01, - .4493018743811285d-01,
     +     0.2729236573866348d-01,0.3547473950556990d-03,
     +     0.1571366799739551d-01,0.4990099219278567d-01,
     +     0.1377915552666770d-01,0.2878206423099872d-02/
c
c***first executable statement d113re
c
c   assign values to w.
c
      do 10 i = 1,13
          do 10 j = 1,5
              w(j,i) = dim3w(i,j)
10    continue
c
c   assign values to g.
c
      do 20 i = 1,3
          do 20 j = 1,13
              g(i,j) = 0
20    continue
      g(1,2) = dim3g(1)
      g(1,3) = dim3g(2)
      g(1,4) = dim3g(3)
      g(1,5) = dim3g(4)
      g(1,6) = dim3g(5)
      g(1,7) = dim3g(6)
      g(2,7) = g(1,7)
      g(1,8) = dim3g(7)
      g(2,8) = g(1,8)
      g(1,9) = dim3g(8)
      g(2,9) = g(1,9)
      g(3,9) = g(1,9)
      g(1,10) = dim3g(9)
      g(2,10) = g(1,10)
      g(3,10) = g(1,10)
      g(1,11) = dim3g(10)
      g(2,11) = g(1,11)
      g(3,11) = g(1,11)
      g(1,12) = dim3g(12)
      g(2,12) = dim3g(11)
      g(3,12) = g(2,12)
      g(1,13) = dim3g(13)
      g(2,13) = g(1,13)
      g(3,13) = dim3g(14)
c
c   assign values to rulpts.
c
      rulpts(1) = 1
      rulpts(2) = 6
      rulpts(3) = 6
      rulpts(4) = 6
      rulpts(5) = 6
      rulpts(6) = 6
      rulpts(7) = 12
      rulpts(8) = 12
      rulpts(9) = 8
      rulpts(10) = 8
      rulpts(11) = 8
      rulpts(12) = 24
      rulpts(13) = 24
c
c   assign values to errcof.
c
      errcof(1) = 4
      errcof(2) = 4.
      errcof(3) = 0.5
      errcof(4) = 3.
      errcof(5) = 0.5
      errcof(6) = 0.25
c
c***end d113re
c
      return
      end
      subroutine d09hre(ndim,wtleng,w,g,errcof,rulpts)
c***begin prologue d09hre
c***keywords basic integration rule, degree 9
c***purpose  to initialize a degree 9 basic rule and null rules.
c***author   alan genz, computer science department, washington
c            state university, pullman, wa 99163-1210 usa
c***last modification 88-05-20
c***description  d09hre initializes a degree 9 integration rule,
c            two degree 7 null rules, one degree 5 null rule and one
c            degree 3 null rule for the hypercube [-1,1]**ndim.
c
c   on entry
c
c   ndim   integer.
c          number of variables.
c   wtleng integer.
c          the number of weights in each of the rules.
c
c   on return
c   w      real array of dimension (5,wtleng).
c          the weights for the basic and null rules.
c          w(1,1),...,w(1,wtleng) are weights for the basic rule.
c          w(i,1),...,w(i,wtleng), for i > 1 are null rule weights.
c   g      real array of dimension (ndim, wtleng).
c          the fully symmetric sum generators for the rules.
c          g(1, j), ..., g(ndim, j) are the are the generators for the
c          points associated with the jth weights.
c   errcof real array of dimension 6.
c          heuristic error coefficients that are used in the
c          error estimation in basrul.
c   rulpts real array of dimension wtleng.
c          a work array.
c
c***references a. genz and a. malik,
c             "an imbedded family of fully symmetric numerical
c              integration rules",
c              siam j numer. anal. 20 (1983), pp. 580-588.
c***routines called-none
c***end prologue d09hre
c
c   global variables
c
      integer ndim,wtleng
      double precision w(5,wtleng),g(ndim,wtleng),errcof(6)
      double precision rulpts(wtleng)
c
c   local variables
c
      double precision ratio,lam0,lam1,lam2,lam3,lamp,twondm
      integer i,j
c
c***first executable statement d09hre
c
c
c     initialize generators, weights and rulpts
c
      do 30 j = 1,wtleng
          do 10 i = 1,ndim
              g(i,j) = 0
10        continue
          do 20 i = 1,5
              w(i,j) = 0
20        continue
          rulpts(j) = 2*ndim
30    continue
      twondm = 2**ndim
      rulpts(wtleng) = twondm
      if (ndim.gt.2) rulpts(8) = (4*ndim* (ndim-1)* (ndim-2))/3
      rulpts(7) = 4*ndim* (ndim-1)
      rulpts(6) = 2*ndim* (ndim-1)
      rulpts(1) = 1
c
c     compute squared generator parameters
c
      lam0 = 0.4707
      lam1 = 4/ (15-5/lam0)
      ratio = (1-lam1/lam0)/27
      lam2 = (5-7*lam1-35*ratio)/ (7-35*lam1/3-35*ratio/lam0)
      ratio = ratio* (1-lam2/lam0)/3
      lam3 = (7-9* (lam2+lam1)+63*lam2*lam1/5-63*ratio)/
     +       (9-63* (lam2+lam1)/5+21*lam2*lam1-63*ratio/lam0)
      lamp = 0.0625
c
c     compute degree 9 rule weights
c
      w(1,wtleng) = 1/ (3*lam0)**4/twondm
      if (ndim.gt.2) w(1,8) = (1-1/ (3*lam0))/ (6*lam1)**3
      w(1,7) = (1-7* (lam0+lam1)/5+7*lam0*lam1/3)/
     +         (84*lam1*lam2* (lam2-lam0)* (lam2-lam1))
      w(1,6) = (1-7* (lam0+lam2)/5+7*lam0*lam2/3)/
     +         (84*lam1*lam1* (lam1-lam0)* (lam1-lam2)) -
     +         w(1,7)*lam2/lam1 - 2* (ndim-2)*w(1,8)
      w(1,4) = (1-9* ((lam0+lam1+lam2)/7- (lam0*lam1+lam0*lam2+
     +         lam1*lam2)/5)-3*lam0*lam1*lam2)/
     +         (18*lam3* (lam3-lam0)* (lam3-lam1)* (lam3-lam2))
      w(1,3) = (1-9* ((lam0+lam1+lam3)/7- (lam0*lam1+lam0*lam3+
     +         lam1*lam3)/5)-3*lam0*lam1*lam3)/
     +         (18*lam2* (lam2-lam0)* (lam2-lam1)* (lam2-lam3)) -
     +         2* (ndim-1)*w(1,7)
      w(1,2) = (1-9* ((lam0+lam2+lam3)/7- (lam0*lam2+lam0*lam3+
     +         lam2*lam3)/5)-3*lam0*lam2*lam3)/
     +         (18*lam1* (lam1-lam0)* (lam1-lam2)* (lam1-lam3)) -
     +         2* (ndim-1)* (w(1,7)+w(1,6)+ (ndim-2)*w(1,8))
c
c     compute weights for 2 degree 7, 1 degree 5 and 1 degree 3 rules
c
      w(2,wtleng) = 1/ (108*lam0**4)/twondm
      if (ndim.gt.2) w(2,8) = (1-27*twondm*w(2,9)*lam0**3)/ (6*lam1)**3
      w(2,7) = (1-5*lam1/3-15*twondm*w(2,wtleng)*lam0**2* (lam0-lam1))/
     +          (60*lam1*lam2* (lam2-lam1))
      w(2,6) = (1-9* (8*lam1*lam2*w(2,7)+twondm*w(2,wtleng)*lam0**2))/
     +         (36*lam1*lam1) - 2*w(2,8)* (ndim-2)
      w(2,4) = (1-7* ((lam1+lam2)/5-lam1*lam2/3+twondm*w(2,
     +         wtleng)*lam0* (lam0-lam1)* (lam0-lam2)))/
     +         (14*lam3* (lam3-lam1)* (lam3-lam2))
      w(2,3) = (1-7* ((lam1+lam3)/5-lam1*lam3/3+twondm*w(2,
     +         wtleng)*lam0* (lam0-lam1)* (lam0-lam3)))/
     +         (14*lam2* (lam2-lam1)* (lam2-lam3)) - 2* (ndim-1)*w(2,7)
      w(2,2) = (1-7* ((lam2+lam3)/5-lam2*lam3/3+twondm*w(2,
     +         wtleng)*lam0* (lam0-lam2)* (lam0-lam3)))/
     +         (14*lam1* (lam1-lam2)* (lam1-lam3)) -
     +         2* (ndim-1)* (w(2,7)+w(2,6)+ (ndim-2)*w(2,8))
      w(3,wtleng) = 5/ (324*lam0**4)/twondm
      if (ndim.gt.2) w(3,8) = (1-27*twondm*w(3,9)*lam0**3)/ (6*lam1)**3
      w(3,7) = (1-5*lam1/3-15*twondm*w(3,wtleng)*lam0**2* (lam0-lam1))/
     +          (60*lam1*lam2* (lam2-lam1))
      w(3,6) = (1-9* (8*lam1*lam2*w(3,7)+twondm*w(3,wtleng)*lam0**2))/
     +         (36*lam1*lam1) - 2*w(3,8)* (ndim-2)
      w(3,5) = (1-7* ((lam1+lam2)/5-lam1*lam2/3+twondm*w(3,
     +         wtleng)*lam0* (lam0-lam1)* (lam0-lam2)))/
     +         (14*lamp* (lamp-lam1)* (lamp-lam2))
      w(3,3) = (1-7* ((lam1+lamp)/5-lam1*lamp/3+twondm*w(3,
     +         wtleng)*lam0* (lam0-lam1)* (lam0-lamp)))/
     +         (14*lam2* (lam2-lam1)* (lam2-lamp)) - 2* (ndim-1)*w(3,7)
      w(3,2) = (1-7* ((lam2+lamp)/5-lam2*lamp/3+twondm*w(3,
     +         wtleng)*lam0* (lam0-lam2)* (lam0-lamp)))/
     +         (14*lam1* (lam1-lam2)* (lam1-lamp)) -
     +         2* (ndim-1)* (w(3,7)+w(3,6)+ (ndim-2)*w(3,8))
      w(4,wtleng) = 2/ (81*lam0**4)/twondm
      if (ndim.gt.2) w(4,8) = (2-27*twondm*w(4,9)*lam0**3)/ (6*lam1)**3
      w(4,7) = (2-15*lam1/9-15*twondm*w(4,wtleng)*lam0* (lam0-lam1))/
     +         (60*lam1*lam2* (lam2-lam1))
      w(4,6) = (1-9* (8*lam1*lam2*w(4,7)+twondm*w(4,wtleng)*lam0**2))/
     +         (36*lam1*lam1) - 2*w(4,8)* (ndim-2)
      w(4,4) = (2-7* ((lam1+lam2)/5-lam1*lam2/3+twondm*w(4,
     +         wtleng)*lam0* (lam0-lam1)* (lam0-lam2)))/
     +         (14*lam3* (lam3-lam1)* (lam3-lam2))
      w(4,3) = (2-7* ((lam1+lam3)/5-lam1*lam3/3+twondm*w(4,
     +         wtleng)*lam0* (lam0-lam1)* (lam0-lam3)))/
     +         (14*lam2* (lam2-lam1)* (lam2-lam3)) - 2* (ndim-1)*w(4,7)
      w(4,2) = (2-7* ((lam2+lam3)/5-lam2*lam3/3+twondm*w(4,
     +         wtleng)*lam0* (lam0-lam2)* (lam0-lam3)))/
     +         (14*lam1* (lam1-lam2)* (lam1-lam3)) -
     +         2* (ndim-1)* (w(4,7)+w(4,6)+ (ndim-2)*w(4,8))
      w(5,2) = 1/ (6*lam1)
c
c     set generator values
c
      lam0 = sqrt(lam0)
      lam1 = sqrt(lam1)
      lam2 = sqrt(lam2)
      lam3 = sqrt(lam3)
      lamp = sqrt(lamp)
      do 40 i = 1,ndim
          g(i,wtleng) = lam0
40    continue
      if (ndim.gt.2) then
          g(1,8) = lam1
          g(2,8) = lam1
          g(3,8) = lam1
      end if
      g(1,7) = lam1
      g(2,7) = lam2
      g(1,6) = lam1
      g(2,6) = lam1
      g(1,5) = lamp
      g(1,4) = lam3
      g(1,3) = lam2
      g(1,2) = lam1
c
c     compute final weight values.
c     the null rule weights are computed from differences between
c     the degree 9 rule weights and lower degree rule weights.
c
      w(1,1) = twondm
      do 70 j = 2,5
          do 50 i = 2,wtleng
              w(j,i) = w(j,i) - w(1,i)
              w(j,1) = w(j,1) - rulpts(i)*w(j,i)
50        continue
70    continue
      do 80 i = 2,wtleng
          w(1,i) = twondm*w(1,i)
          w(1,1) = w(1,1) - rulpts(i)*w(1,i)
80    continue
c
c     set error coefficients
c
      errcof(1) = 5
      errcof(2) = 5
      errcof(3) = 1.
      errcof(4) = 5
      errcof(5) = 0.5
      errcof(6) = 0.25
c
c***end d09hre
c
      end
      subroutine d07hre(ndim,wtleng,w,g,errcof,rulpts)
c***begin prologue d07hre
c***keywords basic integration rule, degree 7
c***purpose  to initialize a degree 7 basic rule, and null rules.
c***author   alan genz, computer science department, washington
c            state university, pullman, wa 99163-1210 usa
c***last modification 88-05-31
c***description  d07hre initializes a degree 7 integration rule,
c            two degree 5 null rules, one degree 3 null rule and one
c            degree 1 null rule for the hypercube [-1,1]**ndim.
c
c   on entry
c
c   ndim   integer.
c          number of variables.
c   wtleng integer.
c          the number of weights in each of the rules.
c          wtleng must be set equal to 6.
c
c   on return
c   w      real array of dimension (5,wtleng).
c          the weights for the basic and null rules.
c          w(1,1),...,w(1,wtleng) are weights for the basic rule.
c          w(i,1),...,w(i,wtleng), for i > 1 are null rule weights.
c   g      real array of dimension (ndim, wtleng).
c          the fully symmetric sum generators for the rules.
c          g(1, j), ..., g(ndim, j) are the are the generators for the
c          points associated with the jth weights.
c   errcof real array of dimension 6.
c          heuristic error coefficients that are used in the
c          error estimation in basrul.
c   rulpts real array of dimension wtleng.
c          a work array.
c
c***references a. genz and a. malik,
c             "an imbedded family of fully symmetric numerical
c              integration rules",
c              siam j numer. anal. 20 (1983), pp. 580-588.
c***routines called-none
c***end prologue d07hre
c
c   global variables
c
      integer ndim,wtleng
      double precision w(5,wtleng),g(ndim,wtleng),errcof(6)
      double precision rulpts(wtleng)
c
c   local variables
c
      double precision ratio,lam0,lam1,lam2,lamp,twondm
      integer i,j
c
c***first executable statement d07hre
c
c
c     initialize generators, weights and rulpts
c
      do 30 j = 1,wtleng
          do 10 i = 1,ndim
              g(i,j) = 0
10        continue
          do 20 i = 1,5
              w(i,j) = 0
20        continue
          rulpts(j) = 2*ndim
30    continue
      twondm = 2**ndim
      rulpts(wtleng) = twondm
      rulpts(wtleng-1) = 2*ndim* (ndim-1)
      rulpts(1) = 1
c
c     compute squared generator parameters
c
      lam0 = 0.4707
      lamp = 0.5625
      lam1 = 4/ (15-5/lam0)
      ratio = (1-lam1/lam0)/27
      lam2 = (5-7*lam1-35*ratio)/ (7-35*lam1/3-35*ratio/lam0)
c
c     compute degree 7 rule weights
c
      w(1,6) = 1/ (3*lam0)**3/twondm
      w(1,5) = (1-5*lam0/3)/ (60* (lam1-lam0)*lam1**2)
      w(1,3) = (1-5*lam2/3-5*twondm*w(1,6)*lam0* (lam0-lam2))/
     +         (10*lam1* (lam1-lam2)) - 2* (ndim-1)*w(1,5)
      w(1,2) = (1-5*lam1/3-5*twondm*w(1,6)*lam0* (lam0-lam1))/
     +         (10*lam2* (lam2-lam1))
c
c     compute weights for 2 degree 5, 1 degree 3 and 1 degree 1 rules
c
      w(2,6) = 1/ (36*lam0**3)/twondm
      w(2,5) = (1-9*twondm*w(2,6)*lam0**2)/ (36*lam1**2)
      w(2,3) = (1-5*lam2/3-5*twondm*w(2,6)*lam0* (lam0-lam2))/
     +         (10*lam1* (lam1-lam2)) - 2* (ndim-1)*w(2,5)
      w(2,2) = (1-5*lam1/3-5*twondm*w(2,6)*lam0* (lam0-lam1))/
     +         (10*lam2* (lam2-lam1))
      w(3,6) = 5/ (108*lam0**3)/twondm
      w(3,5) = (1-9*twondm*w(3,6)*lam0**2)/ (36*lam1**2)
      w(3,3) = (1-5*lamp/3-5*twondm*w(3,6)*lam0* (lam0-lamp))/
     +         (10*lam1* (lam1-lamp)) - 2* (ndim-1)*w(3,5)
      w(3,4) = (1-5*lam1/3-5*twondm*w(3,6)*lam0* (lam0-lam1))/
     +         (10*lamp* (lamp-lam1))
      w(4,6) = 1/ (54*lam0**3)/twondm
      w(4,5) = (1-18*twondm*w(4,6)*lam0**2)/ (72*lam1**2)
      w(4,3) = (1-10*lam2/3-10*twondm*w(4,6)*lam0* (lam0-lam2))/
     +         (20*lam1* (lam1-lam2)) - 2* (ndim-1)*w(4,5)
      w(4,2) = (1-10*lam1/3-10*twondm*w(4,6)*lam0* (lam0-lam1))/
     +         (20*lam2* (lam2-lam1))
c
c     set generator values
c
      lam0 = sqrt(lam0)
      lam1 = sqrt(lam1)
      lam2 = sqrt(lam2)
      lamp = sqrt(lamp)
      do 40 i = 1,ndim
          g(i,wtleng) = lam0
40    continue
      g(1,wtleng-1) = lam1
      g(2,wtleng-1) = lam1
      g(1,wtleng-4) = lam2
      g(1,wtleng-3) = lam1
      g(1,wtleng-2) = lamp
c
c     compute final weight values.
c     the null rule weights are computed from differences between
c     the degree 7 rule weights and lower degree rule weights.
c
      w(1,1) = twondm
      do 70 j = 2,5
          do 50 i = 2,wtleng
              w(j,i) = w(j,i) - w(1,i)
              w(j,1) = w(j,1) - rulpts(i)*w(j,i)
50        continue
70    continue
      do 80 i = 2,wtleng
          w(1,i) = twondm*w(1,i)
          w(1,1) = w(1,1) - rulpts(i)*w(1,i)
80    continue
c
c     set error coefficients
c
      errcof(1) = 5
      errcof(2) = 5
      errcof(3) = 1
      errcof(4) = 5
      errcof(5) = 0.5
      errcof(6) = 0.25
c
c***end d07hre
c
      end
      subroutine drlhre(ndim,center,hwidth,wtleng,g,w,errcof,numfun,
     +                  funsub,scales,norms,x,null,basval,rgnerr,direct)
c***begin prologue drlhre
c***keywords basic numerical integration rule
c***purpose  to compute basic integration rule values.
c***author   alan genz, computer science department, washington
c            state university, pullman, wa 99163-1210 usa
c***last modification 90-02-06
c***description drlhre computes basic integration rule values for a
c            vector of integrands over a hyper-rectangular region.
c            these are estimates for the integrals. drlhre also computes
c            estimates for the errors and determines the coordinate axis
c            where the fourth difference for the integrands is largest.
c
c   on entry
c
c   ndim   integer.
c          number of variables.
c   center real array of dimension ndim.
c          the coordinates for the center of the region.
c   hwidth real array of dimension ndim.
c          hwidth(i) is half of the width of dimension i of the region.
c   wtleng integer.
c          the number of weights in the basic integration rule.
c   g      real array of dimension (ndim,wtleng).
c          the fully symmetric sum generators for the rules.
c          g(1,j), ..., g(ndim,j) are the are the generators for the
c          points associated with the jth weights.
c   w      real array of dimension (5,wtleng).
c          the weights for the basic and null rules.
c          w(1,1),...,w(1,wtleng) are weights for the basic rule.
c          w(i,1),...,w(i,wtleng), for i > 1 are null rule weights.
c   errcof real array of dimension 6.
c          the error coefficients for the rules.
c          it is assumed that the error is computed using:
c           if (n1*errcof(1) < n2 and n2*errcof(2) < n3)
c             then error = errcof(3)*n1
c             else error = errcof(4)*max(n1, n2, n3)
c           error = error + ep*(errcof(5)*error/(es+error)+errcof(6))
c          where n1-n4 are the null rules, ep is the error
c          for the parent
c          subregion and es is the error for the sibling subregion.
c   numfun integer.
c          number of components for the vector integrand.
c   funsub externally declared subroutine.
c          for computing the components of the integrand at a point x.
c          it must have parameters (ndim,x,numfun,funvls).
c           input parameters:
c            x      real array of dimension ndim.
c                   defines the evaluation point.
c            ndim   integer.
c                   number of variables for the integrand.
c            numfun integer.
c                   number of components for the vector integrand.
c           output parameters:
c            funvls real array of dimension numfun.
c                   the components of the integrand at the point x.
c   scales real array of dimension (3,wtleng).
c          scaling factors used to construct new null rules based
c          on a linear combination of two successive null rules
c          in the sequence of null rules.
c   norms  real array of dimension (3,wtleng).
c          2**ndim/(1-norm of the null rule constructed by each of the
c          scaling factors.)
c   x      real array of dimension ndim.
c          a work array.
c   null   real array of dimension (numfun, 8)
c          a work array.
c
c   on return
c
c   basval real array of dimension numfun.
c          the values for the basic rule for each component
c          of the integrand.
c   rgnerr real array of dimension numfun.
c          the error estimates for each component of the integrand.
c   direct real.
c          the coordinate axis where the fourth difference of the
c          integrand values is largest.
c
c***references
c   a.c.genz and a.a.malik, an adaptive algorithm for numerical
c   integration over an n-dimensional rectangular region,
c   j.comp.appl.math., 6:295-302, 1980.
c
c   t.o.espelid, integration rules, null rules and error
c   estimation, reports in informatics 33, dept. of informatics,
c   univ. of bergen, 1988.
c
c***routines called: dfshre, funsub
c
c***end prologue drlhre
c
c   global variables.
c
      external funsub
      integer wtleng,numfun,ndim
      double precision center(ndim),x(ndim),hwidth(ndim),basval(numfun),
     +                 rgnerr(numfun),null(numfun,8),w(5,wtleng),
     +                 g(ndim,wtleng),errcof(6),direct,scales(3,wtleng),
     +                 norms(3,wtleng)
c
c   local variables.
c
      double precision rgnvol,difsum,difmax,frthdf
      integer i,j,k,divaxn
      double precision search,ratio
c
c***first executable statement drlhre
c
c
c       compute volume of subregion, initialize divaxn and rule sums;
c       compute fourth differences and new divaxn (rgnerr is used
c       for a work array here). the integrand values used for the
c       fourth divided differences are accumulated in rule arrays.
c
      rgnvol = 1
      divaxn = 1
      do 10 i = 1,ndim
          rgnvol = rgnvol*hwidth(i)
          x(i) = center(i)
          if (hwidth(i).gt.hwidth(divaxn)) divaxn = i
10    continue
      call funsub(ndim,x,numfun,rgnerr)
      do 30 j = 1,numfun
          basval(j) = w(1,1)*rgnerr(j)
          do 20 k = 1,4
              null(j,k) = w(k+1,1)*rgnerr(j)
20        continue
30    continue
      difmax = 0
      ratio = (g(1,3)/g(1,2))**2
      do 60 i = 1,ndim
          x(i) = center(i) - hwidth(i)*g(1,2)
          call funsub(ndim,x,numfun,null(1,5))
          x(i) = center(i) + hwidth(i)*g(1,2)
          call funsub(ndim,x,numfun,null(1,6))
          x(i) = center(i) - hwidth(i)*g(1,3)
          call funsub(ndim,x,numfun,null(1,7))
          x(i) = center(i) + hwidth(i)*g(1,3)
          call funsub(ndim,x,numfun,null(1,8))
          x(i) = center(i)
          difsum = 0
          do 50 j = 1,numfun
              frthdf = 2* (1-ratio)*rgnerr(j) - (null(j,7)+null(j,8)) +
     +                 ratio* (null(j,5)+null(j,6))
c
c           ignore differences below roundoff
c
              if (rgnerr(j)+frthdf/4.ne.rgnerr(j)) difsum = difsum +
     +            abs(frthdf)
              do 40 k = 1,4
                  null(j,k) = null(j,k) + w(k+1,2)*
     +                        (null(j,5)+null(j,6)) +
     +                        w(k+1,3)* (null(j,7)+null(j,8))
40            continue
              basval(j) = basval(j) + w(1,2)* (null(j,5)+null(j,6)) +
     +                    w(1,3)* (null(j,7)+null(j,8))
50        continue
          if (difsum.gt.difmax) then
              difmax = difsum
              divaxn = i
          end if
60    continue
      direct = divaxn
c
c    finish computing the rule values.
c
      do 90 i = 4,wtleng
          call dfshre(ndim,center,hwidth,x,g(1,i),numfun,funsub,rgnerr,
     +                null(1,5))
          do 80 j = 1,numfun
              basval(j) = basval(j) + w(1,i)*rgnerr(j)
              do 70 k = 1,4
                  null(j,k) = null(j,k) + w(k+1,i)*rgnerr(j)
70            continue
80        continue
90    continue
c
c    compute errors.
c
      do 130 j = 1,numfun
c
c    we search for the null rule, in the linear space spanned by two
c    successive null rules in our sequence, which gives the greatest
c    error estimate among all normalized (1-norm) null rules in this
c    space.
c
          do 110 i = 1,3
              search = 0
              do 100 k = 1,wtleng
                  search = max(search,abs(null(j,i+1)+scales(i,
     +                     k)*null(j,i))*norms(i,k))
100           continue
              null(j,i) = search
110       continue
          if (errcof(1)*null(j,1).le.null(j,2) .and.
     +        errcof(2)*null(j,2).le.null(j,3)) then
              rgnerr(j) = errcof(3)*null(j,1)
          else
              rgnerr(j) = errcof(4)*max(null(j,1),null(j,2),null(j,3))
          end if
          rgnerr(j) = rgnvol*rgnerr(j)
          basval(j) = rgnvol*basval(j)
130   continue
c
c***end drlhre
c
      end
      subroutine dfshre(ndim,center,hwidth,x,g,numfun,funsub,fulsms,
     +                  funvls)
c***begin prologue dfshre
c***keywords fully symmetric sum
c***purpose  to compute fully symmetric basic rule sums
c***author   alan genz, computer science department, washington
c            state university, pullman, wa 99163-1210 usa
c***last modification 88-04-08
c***description dfshre computes a fully symmetric sum for a vector
c            of integrand values over a hyper-rectangular region.
c            the sum is fully symmetric with respect to the center of
c            the region and is taken over all sign changes and
c            permutations of the generators for the sum.
c
c   on entry
c
c   ndim   integer.
c          number of variables.
c   center real array of dimension ndim.
c          the coordinates for the center of the region.
c   hwidth real array of dimension ndim.
c          hwidth(i) is half of the width of dimension i of the region.
c   x      real array of dimension ndim.
c          a work array.
c   g      real array of dimension ndim.
c          the generators for the fully symmetric sum. these must be
c          non-negative and non-increasing.
c   numfun integer.
c          number of components for the vector integrand.
c   funsub externally declared subroutine.
c          for computing the components of the integrand at a point x.
c          it must have parameters (ndim, x, numfun, funvls).
c           input parameters:
c            x      real array of dimension ndim.
c                   defines the evaluation point.
c            ndim   integer.
c                   number of variables for the integrand.
c            numfun integer.
c                   number of components for the vector integrand.
c           output parameters:
c            funvls real array of dimension numfun.
c                   the components of the integrand at the point x.
c   on return
c
c   fulsms real array of dimension numfun.
c          the values for the fully symmetric sums for each component
c          of the integrand.
c   funvls real array of dimension numfun.
c          a work array.
c
c***routines called: funsub
c
c***end prologue dfshre
c
c   global variables.
c
      external funsub
      integer ndim,numfun
      double precision center(ndim),hwidth(ndim),x(ndim),g(ndim),
     +                 fulsms(numfun),funvls(numfun)
c
c   local variables.
c
      integer ixchng,lxchng,i,j,l
      double precision gl,gi
c
c***first executable statement dfshre
c
      do 10 j = 1,numfun
          fulsms(j) = 0
10    continue
c
c     compute centrally symmetric sum for permutation of g
c
20    do 30 i = 1,ndim
          x(i) = center(i) + g(i)*hwidth(i)
30    continue
40    call funsub(ndim,x,numfun,funvls)
      do 50 j = 1,numfun
          fulsms(j) = fulsms(j) + funvls(j)
50    continue
      do 60 i = 1,ndim
          g(i) = - g(i)
          x(i) = center(i) + g(i)*hwidth(i)
          if (g(i).lt.0) go to 40
60    continue
c
c       find next distinct permutation of g and loop back for next sum.
c       permutations are generated in reverse lexicographic order.
c
      do 80 i = 2,ndim
          if (g(i-1).gt.g(i)) then
              gi = g(i)
              ixchng = i - 1
              do 70 l = 1, (i-1)/2
                  gl = g(l)
                  g(l) = g(i-l)
                  g(i-l) = gl
                  if (gl.le.gi) ixchng = ixchng - 1
                  if (g(l).gt.gi) lxchng = l
70            continue
              if (g(ixchng).le.gi) ixchng = lxchng
              g(i) = g(ixchng)
              g(ixchng) = gi
              go to 20
          end if
80    continue
c
c     restore original order to generators
c
      do 90 i = 1,ndim/2
          gi = g(i)
          g(i) = g(ndim-i+1)
          g(ndim-i+1) = gi
90    continue
c
c***end dfshre
c
      end
      subroutine dtrhre(dvflag,ndim,numfun,sbrgns,values,errors,centrs,
     +                  hwidts,greate,error,value,center,hwidth,dir)
c***begin prologue dtrhre
c***purpose dtrhre maintains a heap of subregions.
c***description dtrhre maintains a heap of subregions.
c            the subregions are ordered according to the size
c            of the greatest error estimates of each subregion(greate).
c
c   parameters
c
c     dvflag integer.
c            if dvflag = 1, we remove the subregion with
c            greatest error from the heap.
c            if dvflag = 2, we insert a new subregion in the heap.
c     ndim   integer.
c            number of variables.
c     numfun integer.
c            number of components of the integral.
c     sbrgns integer.
c            number of subregions in the heap.
c     values real array of dimension (numfun,sbrgns).
c            used to store estimated values of the integrals
c            over the subregions.
c     errors real array of dimension (numfun,sbrgns).
c            used to store the corresponding estimated errors.
c     centrs real array of dimension (ndim,sbrgns).
c            used to store the center limits of the stored subregions.
c     hwidts real array of dimension (ndim,sbrgns).
c            used to store the hwidth limits of the stored subregions.
c     greate real array of dimension sbrgns.
c            used to store the greatest estimated errors in
c            all subregions.
c     error  real array of dimension numfun.
c            used as intermediate storage for the error of a subregion.
c     value  real array of dimension numfun.
c            used as intermediate storage for the estimate
c            of the integral over a subregion.
c     center real array of dimension ndim.
c            used as intermediate storage for the center of
c            the subregion.
c     hwidth real array of dimension ndim.
c            used as intermediate storage for the half width of
c            the subregion.
c     dir    integer array of dimension sbrgns.
c            dir is used to store the directions for
c            further subdivision.
c
c***routines called-none
c***end prologue dtrhre
c
c   global variables.
c
      integer dvflag,ndim,numfun,sbrgns
      double precision values(numfun,*),errors(numfun,*)
      double precision centrs(ndim,*)
      double precision hwidts(ndim,*)
      double precision greate(*)
      double precision error(numfun),value(numfun)
      double precision center(ndim),hwidth(ndim)
      double precision dir(*)
c
c   local variables.
c
c   great  is used as intermediate storage for the greatest error of a
c          subregion.
c   direct is used as intermediate storage for the direction of further
c          subdivision.
c   subrgn position of child/parent subregion in the heap.
c   subtmp position of parent/child subregion in the heap.
c
      integer j,subrgn,subtmp
      double precision great,direct
c
c***first executable statement dtrhre
c
c   save values to be stored in their correct place in the heap.
c
      great = greate(sbrgns)
      direct = dir(sbrgns)
      do 5 j = 1,numfun
          error(j) = errors(j,sbrgns)
          value(j) = values(j,sbrgns)
5     continue
      do 10 j = 1,ndim
          center(j) = centrs(j,sbrgns)
          hwidth(j) = hwidts(j,sbrgns)
10    continue
c
c    if dvflag = 1, we will remove the region
c    with greatest estimated error from the heap.
c
      if (dvflag.eq.1) then
          sbrgns = sbrgns - 1
          subrgn = 1
20        subtmp = 2*subrgn
          if (subtmp.le.sbrgns) then
              if (subtmp.ne.sbrgns) then
c
c   find max. of left and right child.
c
                  if (greate(subtmp).lt.greate(subtmp+1)) then
                      subtmp = subtmp + 1
                  end if
              end if
c
c   compare max.child with parent.
c   if parent is max., then done.
c
              if (great.lt.greate(subtmp)) then
c
c   move the values at position subtmp up the heap.
c
                  greate(subrgn) = greate(subtmp)
                  do 25 j = 1,numfun
                      errors(j,subrgn) = errors(j,subtmp)
                      values(j,subrgn) = values(j,subtmp)
25                continue
                  dir(subrgn) = dir(subtmp)
                  do 30 j = 1,ndim
                      centrs(j,subrgn) = centrs(j,subtmp)
                      hwidts(j,subrgn) = hwidts(j,subtmp)
30                continue
                  subrgn = subtmp
                  go to 20
              end if
          end if
      else if (dvflag.eq.2) then
c
c   if dvflag = 2, then insert new region in the heap.
c
          subrgn = sbrgns
40        subtmp = subrgn/2
          if (subtmp.ge.1) then
c
c   compare max.child with parent.
c   if parent is max, then done.
c
              if (great.gt.greate(subtmp)) then
c
c   move the values at position subtmp down the heap.
c
                  greate(subrgn) = greate(subtmp)
                  do 45 j = 1,numfun
                      errors(j,subrgn) = errors(j,subtmp)
                      values(j,subrgn) = values(j,subtmp)
45                continue
                  dir(subrgn) = dir(subtmp)
                  do 50 j = 1,ndim
                      centrs(j,subrgn) = centrs(j,subtmp)
                      hwidts(j,subrgn) = hwidts(j,subtmp)
50                continue
                  subrgn = subtmp
                  go to 40
              end if
          end if
      end if
c
c    insert the saved values in their correct places.
c
      if (sbrgns.gt.0) then
          greate(subrgn) = great
          do 55 j = 1,numfun
              errors(j,subrgn) = error(j)
              values(j,subrgn) = value(j)
55        continue
          dir(subrgn) = direct
          do 60 j = 1,ndim
              centrs(j,subrgn) = center(j)
              hwidts(j,subrgn) = hwidth(j)
60        continue
      end if
c
c***end dtrhre
c
      return
      end
