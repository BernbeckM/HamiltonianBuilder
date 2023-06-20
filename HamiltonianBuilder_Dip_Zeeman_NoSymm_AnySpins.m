%This code builds and diagonalizes the Hamiltonian describing two spins
%interacting via the dipolar interaction, then the Zeeman as successive
%perturbations.

%You will need these variables before you start:

%g1 and g2 - 3x3 diagonal g tensors for your two spins S1 and S2

%g1rottomag and g2rottomag - matrix describing a rotation of your g from
%   crystallographic space into magnetic space. The reason I have it set up
%   this way and not as its transpose (the opposite rotation) is because
%   that's how it's most easily copied from SINGLE_ANISO

%CrystBas1 - how the code is written by default is to use the
%   crystallographic vectors describing the anisotropy axes of your ground
%   multiplet (called a dipole doublet) as the main magnetic axes for which
%   the entire Hamiltonian is described (aka where Bx, By, and Bz will align).
%   Change this definition in the code if you wish. (CrystBas = XXXX)

%First I run another script I just copy/pasted here twice that builds your
%spin operators. It just follows the algorithm from %https://wikipedia.org/wiki/Spin_(physics)

prompt = "What is your first spin? (What is 2 * S1?): ";
spins = input(prompt); %This tells the code the dimensions of your spin operator S1.
Degen1 = spins + 1;

%Create Operators

Sx1 = []; Sy1 = []; Sz1 = [];

%loop for rows

for a = [1:(spins+1)]
    
    %loop for columns

    for b = [1:(spins+1)]
        
        %Kronecker Delta a,b
        if a == b
           
            Sx1(a, b) = 0;
            Sy1(a, b) = 0;
          
            Sz1(a, b) = (spins/2 + 1 - a);
                  
        else

            Sz1(a, b) = 0;

            %handles Kronecker delta for a,b+1 and a+1,b case

            if abs(a - b) == 1

                Sx1(a, b) = sqrt(((spins/2 + 1) * (a + b - 1)) - (a * b))/2;
                Sy1(a, b) = (1i * (a - b) * sqrt(((spins/2 + 1) * (a + b - 1)) - (a * b)))/2;
                
            else

                Sx1(a, b) = 0;
                Sy1(a, b) = 0;
               
            end

        end

    end

end

prompt = "What is your second spin? (What is 2 * S2?): ";
spins = input(prompt); %This tells the code the dimensions of your spin operator S2.
Degen2 = spins + 1;

%Create Operators

Sx2 = []; Sy2 = []; Sz2 = [];

%loop for rows

for a = [1:(spins+1)]
    
    %loop for columns

    for b = [1:(spins+1)]
        
        %Kronecker Delta a,b
        if a == b
           
            Sx2(a, b) = 0;
            Sy2(a, b) = 0;
          
            Sz2(a, b) = (spins/2 + 1 - a);
                  
        else

            Sz2(a, b) = 0;

            %handles Kronecker delta for a,b+1 and a+1,b case

            if abs(a - b) == 1

                Sx2(a, b) = sqrt(((spins/2 + 1) * (a + b - 1)) - (a * b))/2;
                Sy2(a, b) = (1i * (a - b) * sqrt(((spins/2 + 1) * (a + b - 1)) - (a * b)))/2;
                
            else

                Sx2(a, b) = 0;
                Sy2(a, b) = 0;
               
            end

        end

    end

end


%This is requesting the set of fields you want to solve the Zeeman
%Hamiltonian for. The Zeeman and total Hamiltonians will be solved and then
%discarded at each field; only the last is saved. If you have a specific
%field you want to run, just run it as [Field] (1x1 array)
%This builds an array Energies so you can plot(Fields, Energies) and create
%a Zeeman diagram.

disp('Please provide the fields you want in Tesla as a row vector.');
prompt = 'B = '; %Space to make it look nice int he command line

Fields = input(prompt); %This input should be a row vector of the type
                        %[Field1, Field2, .... , FieldLast]. The form
                        %[Field1 : FieldSpacing : FieldLast] also works

Energies = []; %Creating arrays to fill later

B_Plot = []; %I don't know why I built this. Possibly as a backup parameter
             %for if Fields didn't work to plot

CrystBas = CrystBas1; %You need to solve this ahead of time. It's the basis
                      %that rotates FROM diagonal magnetic space TO
                      %cartesian space. It should be an orthonormal basis.
                      %Its transpose will be used. I made it of this form
                      %because that's how I solved for it originally.

g1_DD1_basis = CrystBas' * g1rottomag' * g1; %This creates g tensors for
g2_DD1_basis = CrystBas' * g2rottomag' * g2; %each center using input 
                                             %diagonal g tensors and
                                             %rotation matrices from
                                             %SINGLE_ANISO. It rotates both
                                             %into crystallographic
                                             %coordinates, then it rotates
                                             %both into the diagonal
                                             %magnetic coordinate basis of
                                             %the coupled state you're
                                             %sovling for

%This next bit separates your g tensors into individual variables. I don't
%know why I did it this way - it'd be more efficient to use these indices
%later. I guess I just wanted to be able to look at them, or it made it
%easier to keep track of them in the equations.

gxx1_DD = g1_DD1_basis(1, 1); gxy1_DD = g1_DD1_basis(1, 2); gxz1_DD = g1_DD1_basis(1, 3);
gyx1_DD = g1_DD1_basis(2, 1); gyy1_DD = g1_DD1_basis(2, 2); gyz1_DD = g1_DD1_basis(2, 3);
gzx1_DD = g1_DD1_basis(3, 1); gzy1_DD = g1_DD1_basis(3, 2); gzz1_DD = g1_DD1_basis(3, 3);

gxx2_DD = g2_DD1_basis(1, 1); gxy2_DD = g2_DD1_basis(1, 2); gxz2_DD = g2_DD1_basis(1, 3);
gyx2_DD = g2_DD1_basis(2, 1); gyy2_DD = g2_DD1_basis(2, 2); gyz2_DD = g2_DD1_basis(2, 3);
gzx2_DD = g2_DD1_basis(3, 1); gzy2_DD = g2_DD1_basis(3, 2); gzz2_DD = g2_DD1_basis(3, 3);


%Here I'm defining constants that will be important later.

r_DD = CrystBas' * r; %Internuclear vector rotates to dipole states' magnetic basis
rx = r_DD(1); ry = r_DD(2); rz = r_DD(3); %Breaking r into components (same rant as above)
r_length = norm(r); %Finds the length. Yeah.
C_dip = 0.43297; %This is all the constant terms of the dipolar coupling equatio combined.
uB = 0.4669; %Bohr magneton in cm-1 / T

%This next section sets up the dipolar coupling Hamiltonian. With what
%we've defined beforehand, this takes your g tensors (rotates into the
%output states' magnetic axes), your spin operators, and the internuclear
%vector to build this Hamiltonian. It's split into two parts to make it
%easier to keep track of (so you can check your work later).

%The first term is (S' * g1) * (g2 * S). When you have two terms of S (like
%S1x * S2y), you take their Kronecker product.

Dipole_FirstTerm = ...
      kron((gxx1_DD * Sx1 + gxy1_DD * Sy1 + gxz1_DD * Sz1), (gxx2_DD * Sx2 + gxy2_DD * Sy2 + gxz2_DD * Sz2)) ...
    + kron((gyx1_DD * Sx1 + gyy1_DD * Sy1 + gyz1_DD * Sz1), (gyx2_DD * Sx2 + gyy2_DD * Sy2 + gyz2_DD * Sz2)) ...
    + kron((gzx1_DD * Sx1 + gzy1_DD * Sy1 + gzz1_DD * Sz1), (gzx2_DD * Sx2 + gzy2_DD * Sy2 + gzz2_DD * Sz2));

%This second term is -3/r^2 * (r * g1 * S) * (r * g2 * S). Same rules for
%Kronecker products as above.

Dipole_SecondTerm = (-3 / r_length^2) *...
     (kron(((rx * gxx1_DD + ry * gyx1_DD + rz * gzx1_DD) * Sx1), ((rx * gxx2_DD + ry * gyx2_DD + rz * gzx2_DD) * Sx2)) ...
    + kron(((rx * gxx1_DD + ry * gyx1_DD + rz * gzx1_DD) * Sx1), ((rx * gxy2_DD + ry * gyy2_DD + rz * gzy2_DD) * Sy2)) ...
    + kron(((rx * gxx1_DD + ry * gyx1_DD + rz * gzx1_DD) * Sx1), ((rx * gxz2_DD + ry * gyz2_DD + rz * gzz2_DD) * Sz2)) ...
    ...
    + kron(((rx * gxy1_DD + ry * gyy1_DD + rz * gzy1_DD) * Sy1), ((rx * gxx2_DD + ry * gyx2_DD + rz * gzx2_DD) * Sx2)) ...
    + kron(((rx * gxy1_DD + ry * gyy1_DD + rz * gzy1_DD) * Sy1), ((rx * gxy2_DD + ry * gyy2_DD + rz * gzy2_DD) * Sy2)) ...
    + kron(((rx * gxy1_DD + ry * gyy1_DD + rz * gzy1_DD) * Sy1), ((rx * gxz2_DD + ry * gyz2_DD + rz * gzz2_DD) * Sz2)) ...
    ...
    + kron(((rx * gxz1_DD + ry * gyz1_DD + rz * gzz1_DD) * Sz1), ((rx * gxx2_DD + ry * gyx2_DD + rz * gzx2_DD) * Sx2)) ...
    + kron(((rx * gxz1_DD + ry * gyz1_DD + rz * gzz1_DD) * Sz1), ((rx * gxy2_DD + ry * gyy2_DD + rz * gzy2_DD) * Sy2)) ...
    + kron(((rx * gxz1_DD + ry * gyz1_DD + rz * gzz1_DD) * Sz1), ((rx * gxz2_DD + ry * gyz2_DD + rz * gzz2_DD) * Sz2)));

%This next part combines the two halves of the dipolar Hamiltonian, then
%computes eigenvalues and eigenvectors (column right eigenvectors).

Hamiltonian_Dip = (C_dip / r_length^3) * (Dipole_FirstTerm + Dipole_SecondTerm);
[DipoleVectors, DipoleValues] = eig(Hamiltonian_Dip);

%This part builds your magnetic transition matrices from the Zeeman
%Hamiltonian. It is the derivative of the Zeeman Hamiltonian with respect
%to a field Bi, which is the field B aligned entirely along axis i.
%(i = X, Y, and Z)

Transition_X = (kron((gxx1_DD * Sx1 + gxy1_DD * Sy1 + gxz1_DD * Sz1), eye(Degen2)) + kron(eye(Degen1), (gxx2_DD * Sx2 + gxy2_DD * Sy2 + gxz2_DD * Sz2)));
Transition_Y = (kron((gyx1_DD * Sx1 + gyy1_DD * Sy1 + gyz1_DD * Sz1), eye(Degen2)) + kron(eye(Degen1), (gyx2_DD * Sx2 + gyy2_DD * Sy2 + gyz2_DD * Sz2)));
Transition_Z = (kron((gzx1_DD * Sx1 + gzy1_DD * Sy1 + gzz1_DD * Sz1), eye(Degen2)) + kron(eye(Degen1), (gzx2_DD * Sx2 + gzy2_DD * Sy2 + gzz2_DD * Sz2)));

%This next finds the average transition matrix using the method described
%in OpenMolcas, namely (T_X + T_Y + T_Z)/3

Transition_Avg = (abs(Transition_X) + abs(Transition_Y) + abs(Transition_X)) ./ 3;

%This next section builds the Zeeman Hamiltonian proper, with the results
%computed over the entire array Fields

for B_tot = Fields 

%Define here how you want to align your field. These axes correspond to the
%cartesian axes described by the diagonal g tensor of your chosen reference
%state (typically the ground state, DD1)

Bx = 0;
By = 0;
Bz = B_tot;

%The Zeeman Hamiltonian is built in steps. The equation for an uncoupled spin is:
%[Bx By Bz] [gi] [S]. This yields a scalar multiple of your chosen 
%pin matrices that we can interact later.

gS1_Zeem = Bx * (gxx1_DD * Sx1 + gxy1_DD * Sy1 + gxz1_DD * Sz1)...
         + By * (gyx1_DD * Sx1 + gyy1_DD * Sy1 + gyz1_DD * Sz1)...
         + Bz * (gzx1_DD * Sx1 + gzy1_DD * Sy1 + gzz1_DD * Sz1);

gS2_Zeem = Bx * (gxx2_DD * Sx2 + gxy2_DD * Sy2 + gxz2_DD * Sz2)...
         + By * (gyx2_DD * Sx2 + gyy2_DD * Sy2 + gyz2_DD * Sz2)...
         + Bz * (gzx2_DD * Sx2 + gzy2_DD * Sy2 + gzz2_DD * Sz2);

%The second step is to find the Kronecker sum of gS1 and gS2. This function
%isn't built into matlab, so we set it up manually here. This is why Degen
%was defined before - to build identity matrices of the correct dimensions.

Hamiltonian_Zeeman = uB * (kron(gS1_Zeem, eye(Degen2)) + kron(eye(Degen1), gS2_Zeem));

%Here, the eigenvalues and right column eigenvectors are found.

[ZeemanVectors, ZeemanValues] = eig(Hamiltonian_Zeeman);

%Finally, we combine the two Hamiltonians and find the eigenvalues and
%eigenvectors of the total Hamiltonian. This should be different than the
%sum of the two solved separately.

Hamiltonian_Total = Hamiltonian_Dip + Hamiltonian_Zeeman;
[TotalVectors, TotalValues] = eig(Hamiltonian_Total);

%Since we're doing a loop, this stores the values so we can plot them
%later. It saves only the eigenvalues though - vectors are disarded. If you
%want a specific set of eigenvalues and eigenvectors, use a 1x1 array.

Energies =  [Energies, diag(TotalValues)]; 
B_Plot = [B_Plot, B_tot * ones(length(diag(TotalValues)),1)]; 
%Fields = B_Plot;
end