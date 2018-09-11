% Figure with the plots Outstanding Debt in Steady State
% Ranges for debt levels
dt = parameters.dt;
% Range for years. Up to 1 year. Up to 2 years etc. 
% Note that they start with 0.5 years, because its July!
years           = [0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5 14.5 17.5 19.5];

calendar_years  = ['2018' '2019' '2020' '2021' '2022' '2023' '2024' '2025' '2026' '2027' '2028'...
                   '2029' '2030' '2031/32' '2033/36' '2037/39'];

debt_data       = [5.7  13.07  8.68 7.1  6.79  5.93   8.20  7.22  6.84  5.38 4.93    2.52   3.01   2.22   1.62    2.0]; % This data comes from the Spanish Treasury, July 2018

for i=1:(size(years,2))
    if i==1
    ranges{i}= (1:years(i)/dt);
    else 
    ranges{i}= (years(i-1)/dt+1:years(i)/dt);
    end
end

% Building Densities
for i=1:size(years,2)
f_i_rss (i) = sum(results.f_rss(ranges{i},end))*dt;
end

%comparable numbers
f_i_rss_pct   = f_i_rss/sum(f_i_rss);

string_calendar_years  = {'2018' '2019' '2020' '2021' '2022' '2023' '2024' '2025' '2026'  '2027' '2028'...
                   '2029'  '2030'  '2031/32'  '2033/36' '2037/39'};
               
bar_data = [debt_data_pct' f_i_rss_pct'];
bar(years,bar_data)
set(gca, 'XTickLabel',string_calendar_years, 'XTick',1:numel(string_calendar_years))

bar([f_i_rss_pct' debt_data_pct'])
set(gca, 'XTickLabel',string_calendar_years, 'XTick',1:numel(string_calendar_years))
xlabel('Maturity','interpreter','LaTex','FontSize',12)
ylabel('$\%$ of GDP','interpreter','LaTex','FontSize',12)
legend({'Model','Data'},'interpreter','LaTex','FontSize',12)


% cd figures; save Maturity_Profile; figure(2); saveas(gcf,'Maturity_Profile','pdf');
