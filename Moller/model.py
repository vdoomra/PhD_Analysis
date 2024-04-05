    model = build_model_two_hidden_layers()
    model.summary()

    EPOCHS = 200
    batch_size = 32

    history = 

    plotter = tfdocs.plots.HistoryPlotter(smoothing_std=2)
    plotter.plot({'Basic' : history}, metric = "mse")
    plt.ylabel('MAE')
    plt.show()

    residuals = predict_in_batches(model, x_test_scaled, y_test, batch_size=32)

    plt.figure(figsize=(10, 6))
    plt.hist(residuals, bins=20, density=True, alpha=0.7)
    plt.title('Histogram of Differences between Actual and Predicted Values')
    plt.xlim([-0.02,0.02])
    plt.xlabel('Difference')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)

    plt.show()