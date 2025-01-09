function err = compareSurrogate(m_s, m_ode45,T)
    
    t = (m_ode45-m_s).^2;
    err = (1/length(T))*trapz(t);
    
end