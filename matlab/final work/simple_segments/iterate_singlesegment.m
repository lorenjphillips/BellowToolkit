function ss_outputs = iterate_singlesegment(max_theta,min_ell,max_ell,step_ell)
    % Get the range and step size for ell from the user
    ell_start = min_ell;%input('Enter the start of the range for ell: ');
    ell_end = max_ell;%input('Enter the end of the range for ell: ');
    ell_step = step_ell;%input('Enter the step size for ell: ');
    
    % Generate the ell values
    ell_values = ell_start:ell_step:ell_end;
    numEll = length(ell_values);
    
    % Preallocate the output 3D matrix
    numSteps = max_theta + 1;
    ss_outputs = zeros(numSteps, 2, numEll);
    
    % Loop through ell values and calculate the ss matrix for each
    for i = 1:numEll
        ell = ell_values(i);
        ss = singlesegment_nooutput(max_theta, ell);
        ss_outputs(:, :, i) = ss;
    end
    
    % Display the results as needed
    % disp('3D Matrix of Angles, Distances, and Ell values:');
    % disp(output_3D_matrix);
    fprintf(2,'Success. ')
end
