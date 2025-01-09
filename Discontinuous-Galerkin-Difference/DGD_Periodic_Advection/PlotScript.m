function PlotScript(phi)

for i=1:length(phi)
    figure
    for j=1:length(phi(i).x)
        plot(phi(i).x{j},phi(i).val{j});
        hold on;
        grid on;
    end
end

% for i=1:length(phi)
%     figure
%     for j=1:length(phi(i).x)
%         plot(phi(i).x{j},phi(i).der{j});
%         hold on;
%         grid on;
%     end
% end

end