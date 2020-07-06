function [Floe] = Weld_Floes(Floe,dhdt,Amax,SUBFLOES)
%%This function takes in two floe field and based upon a specified
%%probability function welds interacting floes together

i = 1; 
floenew = [];

%Set probability function
ramp = @(frac) heaviside(frac)*frac;
SimpMin = @(A) 15*log10(A);

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
                        floe = FuseFloes(Floe(i),Floe(j),SUBFLOES);
                        A(i) = Floe(i).area;
                        Floe(j).alive  = 0;
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

