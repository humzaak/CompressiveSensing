for i = 1:50
    error(i) = test_OMP_and_CoSaMP(i);
end
plot(error);
xlabel('SNR in dB');
ylabel('norm(xr-x)/norm(x)');