artlessV2 <-
function(dat, z, x=NULL, pr=NULL, xm = NULL, near = NULL, fine = NULL,
                    ncontrols = 1, rnd = 2, solver = "rlemon"){
  #  Call alittleArt using its default penalties
  alittleArt(dat,z,x=x,pr=pr,xm=xm,near=near,fine=fine,ncontrols=ncontrols,
             rnd=rnd,solver=solver,min.penalty=c(1,1,0.05),
             pr.penalty=c(1,2,25,250))
}
