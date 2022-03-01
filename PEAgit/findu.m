function val = findu(u0,mp,k,n,z,theta)

val = mp.alpha * theta * u0^(mp.alpha-1)*k^(mp.alpha-1)*n^(1-mp.alpha) - z*mp.phi1 - mp.phi2*(u0-1);