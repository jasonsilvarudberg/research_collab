% Create adjusted gray matter and CSF probability maps to address voxel averaging effects in diffusion weighted imaging. Based on equations:
% GM signal:
% signal_gm = rho_gm * exp(-TE / T2_gm) * (1 - exp(-TR / T1_gm))% 
% 
% 
% CSF signal:
% signal_csf = rho_csf * exp(-TE / T2_csf) * (1 - exp(-TR / T1_csf))
% 
% 
% Weighted GM:
% weighted_gm = lambda_gm * signal_gm
% 
% 
% Weighted CSF:
% weighted_csf = lambda_csf * signal_csf
% 
% 
% Adjusted GM probability:
% adjusted_gm = weighted_gm / (weighted_gm + weighted_csf)
% 
% 
% Adjusted CSF probability:
% adjusted_csf = weighted_csf / (weighted_gm + weighted_csf)
%
% For each subject, this script:
%   1. Reads the registered GM and CSF probability maps.
%   2. Calculates the expected GM and CSF MRI signal.
%   3. Weights the tissue probability maps by the expected signal.
%   4. Normalizes the weighted maps voxel by voxel.
%   5. Saves adjusted GM and CSF NIfTI files.
%
% Variable definitions:
%   lambda_gm  = GM probability map in diffusion space
%   lambda_csf = CSF probability map in diffusion space
%   rho_gm     = GM proton density
%   rho_csf    = CSF proton density
%   TE         = echo time
%   TR         = repetition time
%   T1_gm      = GM T1 relaxation time
%   T2_gm      = GM T2 relaxation time
%   T1_csf     = CSF T1 relaxation time
%   T2_csf     = CSF T2 relaxation time


% -------------------------------------------------------------------------
% MRI and tissue parameters
% -------------------------------------------------------------------------
% Replace these values with the values used for your acquisition and model.

rho_csf = j;
TE      = k;
TR      = l;
T1_csf  = m;
T2_csf  = n;
T1_gm   = o;
T2_gm   = p;
rho_gm  = q;


% -------------------------------------------------------------------------
% Subject IDs
% -------------------------------------------------------------------------

subject_ids = {
    'js000'
    ...
};


% -------------------------------------------------------------------------
% Input and output folders
% -------------------------------------------------------------------------
% Use pwd to process files in the current working directory.

input_folder  = pwd;
output_folder = pwd;

if ~isfolder(output_folder)
    mkdir(output_folder);
end


% -------------------------------------------------------------------------
% Calculate expected tissue signal
% -------------------------------------------------------------------------
% Signal = rho * exp(-TE/T2) * (1 - exp(-TR/T1))

signal_gm = rho_gm .* exp(-TE ./ T2_gm) .* ...
    (1 - exp(-TR ./ T1_gm));

signal_csf = rho_csf .* exp(-TE ./ T2_csf) .* ...
    (1 - exp(-TR ./ T1_csf));

fprintf('Expected GM signal:  %.6f\n', signal_gm);
fprintf('Expected CSF signal: %.6f\n\n', signal_csf);


% -------------------------------------------------------------------------
% Process each subject
% -------------------------------------------------------------------------

for subject_index = 1:numel(subject_ids)

    subject_id = subject_ids{subject_index};

    fprintf(
        '[%d/%d] Processing %s...\n', ...
        subject_index, ...
        numel(subject_ids), ...
        subject_id ...
    );

    try
        % Construct input file names.
        gm_base = fullfile(
            input_folder, ...
            [subject_id '_gm_reg']
        );

        csf_base = fullfile(
            input_folder, ...
            [subject_id '_csf_reg']
        );

        % Find either the .nii or .nii.gz version.
        gmprob  = find_nifti_file(gm_base);
        csfprob = find_nifti_file(csf_base);

        % Construct output file names.
        appgm = fullfile(
            output_folder, ...
            [subject_id '_appgm.nii']
        );

        appcsf = fullfile(
            output_folder, ...
            [subject_id '_appcsf.nii']
        );

        % Read the probability maps as double-precision arrays.
        lambda_gm  = double(niftiread(gmprob));
        lambda_csf = double(niftiread(csfprob));

        % Confirm that both images have the same dimensions.
        if ~isequal(size(lambda_gm), size(lambda_csf))
            warning(
                ['Skipping %s because the GM and CSF maps have ' ...
                 'different dimensions.'], ...
                subject_id ...
            );
            continue;
        end

        % Read image metadata.
        gm_info  = niftiinfo(gmprob);
        csf_info = niftiinfo(csfprob);

        % Confirm that the voxel dimensions match.
        if ~isequal(gm_info.PixelDimensions, csf_info.PixelDimensions)
            warning(
                ['Skipping %s because the GM and CSF maps have ' ...
                 'different voxel dimensions.'], ...
                subject_id ...
            );
            continue;
        end

        % Weight each tissue probability map by its expected MRI signal.
        weighted_gm  = lambda_gm .* signal_gm;
        weighted_csf = lambda_csf .* signal_csf;

        % Calculate the voxelwise normalization denominator.
        denominator = weighted_gm + weighted_csf;

        % Initialize adjusted maps as zero.
        adjusted_gm  = zeros(size(lambda_gm), 'double');
        adjusted_csf = zeros(size(lambda_csf), 'double');

        % Normalize only where the denominator is finite and greater than 0.
        valid_voxels = ...
            isfinite(denominator) & ...
            denominator > 0;

        adjusted_gm(valid_voxels) = ...
            weighted_gm(valid_voxels) ./ ...
            denominator(valid_voxels);

        adjusted_csf(valid_voxels) = ...
            weighted_csf(valid_voxels) ./ ...
            denominator(valid_voxels);

        % Replace any remaining nonfinite values with zero.
        adjusted_gm(~isfinite(adjusted_gm)) = 0;
        adjusted_csf(~isfinite(adjusted_csf)) = 0;

        % Save as single precision to reduce output file size.
        adjusted_gm  = single(adjusted_gm);
        adjusted_csf = single(adjusted_csf);

        % Update metadata for single-precision output.
        output_info = gm_info;
        output_info.Datatype = 'single';
        output_info.BitsPerPixel = 32;
        output_info.Filename = '';

        % Write uncompressed NIfTI files.
        niftiwrite(
            adjusted_gm, ...
            appgm, ...
            output_info, ...
            'Compressed', ...
            false ...
        );

        niftiwrite(
            adjusted_csf, ...
            appcsf, ...
            output_info, ...
            'Compressed', ...
            false ...
        );

        fprintf('    Saved: %s\n', appgm);
        fprintf('    Saved: %s\n\n', appcsf);

    catch ME
        warning(
            'Could not process %s: %s', ...
            subject_id, ...
            ME.message ...
        );
    end
end

fprintf('Processing complete.\n');


% -------------------------------------------------------------------------
% Local function
% -------------------------------------------------------------------------

function file_path = find_nifti_file(file_base)
% Find a NIfTI file stored as either .nii or .nii.gz.
%
% Input:
%   file_base: Full path without a file extension
%
% Output:
%   file_path: Full path to the existing NIfTI file

    nii_file = [file_base '.nii'];
    nii_gz_file = [file_base '.nii.gz'];

    if isfile(nii_file)
        file_path = nii_file;

    elseif isfile(nii_gz_file)
        file_path = nii_gz_file;

    else
        error(
            'Neither "%s" nor "%s" was found.', ...
            nii_file, ...
            nii_gz_file ...
        );
    end
end