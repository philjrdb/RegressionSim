function nullZ_data = nullZ(data)

nullZ_data = data./sqrt(mean(data.^2));