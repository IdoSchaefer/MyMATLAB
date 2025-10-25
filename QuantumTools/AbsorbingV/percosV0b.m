function performance = percosV0b(Vk, xdomain, kdomain, Nk, penalty)
% The function computes the performance parameter for a cosine series
% absorbing potential, with 0 boundary conditions in both edges.
    sz = size(Vk);
    if sz(1) == 1
        Nfree_k = sz(2)/2;
        Vk = Vk.';
    else
        Nfree_k = sz(1)/2;
    end
    % The cosine term coefficients, multiplied by sqrt(Nk/2), for the real
    % and imaginary parts:
    Vk_re_fix0 = Vk(1:Nfree_k);
    Vk_re_fix0(1) = Vk_re_fix0(1)/2;
    Vk_im_fix0 = Vk((Nfree_k + 1):2*Nfree_k);
    Vk_im_fix0(1) = Vk_im_fix0(1)/2;
    % The Nk and Nk-1 terms, subject to the boundary constraints:
    VkNre = -2*sum(Vk_re_fix0(Nfree_k:-2:1));
    VkNminus1re = -sum(Vk_re_fix0((Nfree_k - 1):-2:1));
    VkNim = -2*sum(Vk_im_fix0(Nfree_k:-2:1));
    VkNminus1im = -sum(Vk_im_fix0((Nfree_k - 1):-2:1));
    Vx = [dctI([Vk(1:(Nfree_k)); VkNminus1re; VkNre]) + 1i*dctI([Vk((Nfree_k + 1):2*Nfree_k); VkNminus1im; VkNim])];
    performance = Vabs_efficiency1(Vx, xdomain, kdomain, Nk) + penalty*sum(Vk.^2);
end