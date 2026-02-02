function R = quat_to_rotm(q)
%QUAT_TO_ROTM Convert quaternion q=[qw;qx;qy;qz] to rotation matrix.
q = q / norm(q);
qw = q(1); qx = q(2); qy = q(3); qz = q(4);
R = [1-2*(qy^2+qz^2), 2*(qx*qy-qz*qw), 2*(qx*qz+qy*qw);
     2*(qx*qy+qz*qw), 1-2*(qx^2+qz^2), 2*(qy*qz-qx*qw);
     2*(qx*qz-qy*qw), 2*(qy*qz+qx*qw), 1-2*(qx^2+qy^2)];
end
