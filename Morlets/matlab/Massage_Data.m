function [] = Massage_Data( file_name )

    fileID = fopen(file_name, 'r');

    A = fscanf(fileID, '%f');

    sample_length = size(A, 1) / 3;

    aPrime = reshape(A, sample_length, 3);

    xMat = aPrime(:,1);
    yMat = aPrime(:,2);
    zMat = aPrime(:,3);
    
    file_name_x = strcat(file_name, 'x');
    file_name_x = strcat(file_name_x, '.txt');
    
    file_name_y = strcat(file_name, 'y');
    file_name_y = strcat(file_name_y, '.txt');
    
    file_name_z = strcat(file_name, 'z');
    file_name_z = strcat(file_name_z, '.txt');

    save(file_name_x, 'xMat', '-ascii');
    save(file_name_y, 'yMat', '-ascii');
    save(file_name_z, 'zMat', '-ascii');