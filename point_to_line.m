% calculates the distance between a dot and a line
function dist = point_to_line(pt, v1, v2)
    a = v1 - v2;
    b = pt - v2;
    line_vec = a; %vector(start, end) # (3.5, 0, -1.5)
    pnt_vec = b;  %vector(start, pnt)  # (1, 0, -1.5)
    line_len = sqrt(sum(line_vec.^2)); % # 3.808
    line_unitvec = line_vec/line_len; % # (0.919, 0.0, -0.394)
    pnt_vec_scaled = pnt_vec/line_len; %  # (0.263, 0.0, -0.393)
    t = dot(line_unitvec, pnt_vec_scaled); % # 0.397

    if t < 0.0
         t = 0.0;
    elseif t > 1.0
         t = 1.0;
    end

    nearest = line_vec* t; %    # (1.388, 0.0, -0.595)
    dist = sqrt(sum((nearest-pnt_vec).^2)); % # 0.985
    nearest = nearest+v2;
end