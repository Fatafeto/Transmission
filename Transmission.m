channels = 100;
E_User = 1/144;     %Erlang per user
Pb = 0.005;     %blocking probability
Area = 450*10^6;
k = 2;      %user density per meter squared
user_sens = db2pow(-120);       %user sensitivity
P_Max_In = 1;
P_Max_Out = 2;
C_I_required = 3;        %minimum carrier to interference ratio
lambda = (3*10^8)/(1.8*10^9);
constant = 1.5*sqrt(3);

counter = 1;

for x = setdiff(1:13,[2,5,6,8,10,11])
    for y = setdiff(1:13,[2,5,6,8,10,11])
        N_in = x;       %inner cell reuse factor
        N_out = y;      %outer cell reuse factor
        for i = 1:100
            Channels_in = i;
            Channels_out = 100 - i;
            
            E_in = findrhob(floor(Channels_in/N_in),Pb);        %Erlang increases as the number of trunks increases
            E_out = findrhob(floor(Channels_out/N_out),Pb);     %number of trunks increases if N decreases or number of channels increases or both
            
            R_in_squared = E_in/(k*E_User*pi);
            R_in = sqrt(R_in_squared);
            
            R_out_squared = (E_out + E_in)/(k*E_User*constant);     %R decreases when Erlang decreases
            R_out = sqrt(R_out_squared);
            
            P_sens_in = P_Max_In*((lambda/(4*pi*R_in))^2);
            P_sens_out = P_Max_Out*((lambda/(4*pi*R_out))^2);
            
            C_I_out = (3*N_out)/6;
            C_I_in = (3*R_out_squared*N_in)/(6*R_in_squared);
            
            if P_sens_in >= user_sens && P_sens_out >= user_sens && C_I_in >= C_I_required && C_I_out >= C_I_required && N_out >= N_in
                N_In(1,counter) = N_in;
                N_Out(1,counter) = N_out;
                Diameter_Out(1,counter) = 2*R_out;
                Diameter_In(1,counter) = 2*R_in;
                Channels_In(1,counter) = Channels_in;
                Channels_Out(1,counter) = Channels_out;
                Number_Cells(1,counter) = ceil(Area/(constant*R_out_squared));      %number of cells increases when radius squared decreases
                counter = counter + 1;
                break;
            end
        end
    end
end

%Required
Number_Cells_Min = min(Number_Cells);
index = find(Number_Cells == Number_Cells_Min);
N_in = N_In(index);
N_out = N_Out(index);
Num_Inner_Channels = Channels_In(index);
Num_Outer_Channels = Channels_Out(index);
Diameter_in = Diameter_In(index);
Diameter_out = Diameter_Out(index);

