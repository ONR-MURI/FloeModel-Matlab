function [Floe] = MeldFloes(Floe,dhdt,Amax,SUBFLOES,Nx,Ny,c2_boundary,PERIODIC)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
FloeOld = Floe;
i = 1; 
floenew = [];
ramp = @(frac) heaviside(frac)*frac;
SimpMin = @(A) 15*log10(A);
while i < length(Floe)
    j = length(Floe);
    if Floe(i).alive
        while j > i
            if Floe(j).alive
                
                if sqrt((Floe(i).Xi-Floe(j).Xi)^2 + (Floe(i).Yi-Floe(j).Yi)^2)<(Floe(i).rmax+Floe(j).rmax) % if floes are potentially overlapping
                    polyout = union(Floe(i).poly,Floe(j).poly);
                    areaPoly = area(polyout);
                    frac = (Floe(i).area+Floe(j).area)/Amax;
                    p = rand(1);
                    if areaPoly/Floe(i).area > 0.99
                        p=0;
                    elseif areaPoly/Floe(j).area > 0.99
                       p=0;
                    end
                    if p <ramp((1-frac)^5*dhdt) && polyout.NumRegions ==1 %areaPoly>0
                        floeold = Floe(i);
                        floe2 = Floe(j);
                        floe = FuseFloes(Floe(i),Floe(j),SUBFLOES);
                        A(i) = Floe(i).area;
                        Floe(j).alive  = 0;
                        floe2old = floe;
                        if length(floe.poly.Vertices) > SimpMin(Floe(i).area)
                            floe = FloeSimplify(floe, 250,SUBFLOES);
                        end
                        for jj = 1:length(floe)
                            if jj == 1
                                Floe(i) = floe(jj);
                            else
                                floenew = [floenew floe(jj)];
                            end
                        end
                        

                    end
                
                end
            end
            j = j-1;
            
        end
        
    end
    i = i+1;
end
Floe = [Floe floenew];
end

