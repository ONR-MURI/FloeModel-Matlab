function [Floe] = Weld_Floes(Floe,Nb,dhdt,Amax)
%%This function takes in two floe field and based upon a specified
%%probability function welds interacting floes together

id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

i = Nb+1; 
floenew = [];

for ii =1:length(Floe)
    Floe(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
end

%Set probability function
hvsd = @(x) [0.5*(x == 0) + (x > 0)];
ramp = @(frac) hvsd(frac)*frac;
SimpMin = @(A) 3*log10(A);

%Loop through all floes to determine which will weld
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

                    %If probability is met then weld them together
                    if p <ramp((1-frac)^5*dhdt) && polyout.NumRegions ==1
                        floe = FuseFloes(Floe(i),Floe(j));
                        if (Floe(i).mass+Floe(j).mass)/floe.mass - 1 > 0.001
                            xx = 1;
                            xx(1) = [1 2];
                        end
                        A(i) = Floe(i).area;
                        Floe(j).alive  = 0;
                        if length(floe.poly.Vertices) > SimpMin(Floe(i).area)
                            floe = FloeSimplify(floe);
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
live = cat(1,Floe.alive);
Floe(live == 0) = [];
Floe=rmfield(Floe,{'poly'});

warning('on',id)
warning('on',id3)

end

