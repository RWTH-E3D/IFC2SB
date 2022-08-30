#include "ViewerIncludes.h"


class DialogTransparency : public QDialog {
Q_OBJECT
public:
    explicit DialogTransparency(QWidget *parent = nullptr, Qt::WindowFlags f = nullptr, bool modal = true);

    ~DialogTransparency() override;

signals:

    void sendTransparencyChanged(int value);
};
