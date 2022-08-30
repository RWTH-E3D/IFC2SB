#ifndef MATERIAL_H
#define MATERIAL_H

#include "ViewerIncludes.h"

class DialogMaterial : public QDialog {
Q_OBJECT
public:
    explicit DialogMaterial(QWidget *parent = nullptr, bool modal = true, Qt::WindowFlags f = nullptr);

    ~DialogMaterial() override;

signals:

    void sendMaterialChanged(int);

public slots:

    void updateButtons(bool isOn);

private:
    QList<QPushButton *> myButtons;
};

#endif
