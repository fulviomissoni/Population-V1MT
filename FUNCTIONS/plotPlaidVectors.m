function plotPlaidVectors(truetheta_vec, theta_g_mat, vpld_vec, vgrat_mat)
    % truetheta_vec : [N x 1] plaid directions (rad)
    % theta_g_mat   : [N x 2] grating normal directions (rad)
    % vpld_vec      : [N x 1] plaid speeds (pixels/frame)
    % vgrat_mat     : [N x 2] grating normal speeds (pixels/frame)
    N = numel(truetheta_vec);
    figure; clf;
    set(gcf, 'color', 'w');
    for i = 1:N
        subplot(N, 1, i); hold on; axis equal;
        set(gca, 'XAxisLocation','origin','YAxisLocation','origin');
        axis([-1 1 -1 1] * max(vpld_vec)*1.3);
        truetheta = truetheta_vec(i);
        vpld = vpld_vec(i);
        theta_g = theta_g_mat(i, :);
        vgrat = vgrat_mat(i, :);
        % Plaid velocity vector
        v_plaid = [vpld * cos(truetheta); vpld * sin(truetheta)];
        quiver(0, 0, v_plaid(1), v_plaid(2), 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        % Grating velocities
        colors = ['b','g'];
        for j = 1:2
            n_j = [cos(theta_g(j)); sin(theta_g(j))];
            vj = vgrat(j) * n_j;
            quiver(0, 0, vj(1), vj(2), 0, colors(j), 'LineWidth',1.5, 'MaxHeadSize',0.5);
            % Also draw the normal direction as dashed
            quiver(0, 0, n_j(1), n_j(2), 0, [colors(j) '--']);
        end
        title(sprintf('Case %d', i));
    end
end