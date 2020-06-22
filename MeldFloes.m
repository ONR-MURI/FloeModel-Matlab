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
%                     polyout = intersect(Floe(i).poly,Floe(j).poly);
                    polyout = union(Floe(i).poly,Floe(j).poly);
                    areaPoly = area(polyout);
%                     if areaPoly==0 && polyout2.NumRegions ==1
%                         x = 1;
%                         x(1) = [1 2];
%                     end
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
                        if Floe(i).poly.NumRegions>1
                            xx = 1;
                            xx(1) = [1 2];
                        end
                        A(i) = Floe(i).area;
                        Floe(j).alive  = 0;
                        floe2old = floe;
                        if length(floe.poly.Vertices) > SimpMin(Floe(i).area)
                            floe = FloeSimplify(floe, 250,SUBFLOES);
                        end
                        for jj = 1:length(floe)
                            if floe(jj).poly.NumRegions > 1
                                xx = 1;
                                xx(1) = [1 2];
                            end
                            if jj == 1
                                Floe(i) = floe(jj);
                            else
                                floenew = [floenew floe(jj)];
                            end
                        end
                        
                        
                        [eularian_data] = calc_eulerian_data2([Floe floenew],Nx,Ny,c2_boundary,PERIODIC);
                        if min(min(eularian_data.c))<0.25 || max(max(eularian_data.c))>1.05
                            xx = 1;
                            xx(1) = [1 2];
                        end
%                         for kk = 1:length(Floe)
%                             polyI = intersect(Floe(i).poly,Floe(kk).poly);
%                             if i == kk
%                                 1;
%                             elseif area(polyI)/Floe(kk).area > 0.99
%                                 xx = 1;
%                                 xx(1) = [1 2];
%                             elseif area(polyI)/Floe(i).area > 0.99
%                                 xx = 1;
%                                 xx(1) = [1 2];
%                             end
%                         end
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

