#ifndef TRANSLATE_H
#define TRANSLATE_H

#include "IESample.h"

#include <Standard_WarningsDisable.hxx>
#include <QObject>
#include <Standard_WarningsRestore.hxx>

#include <AIS_InteractiveContext.hxx>
#include <TopTools_HSequenceOfShape.hxx>

class TranslateDlg;

class IESAMPLE_EXPORT Translate : public QObject {
Q_OBJECT

public:
    enum {
        FormatBREP, FormatIGES, FormatSTEP, FormatVRML, FormatSTL
    };

    explicit Translate(QObject *);

    ~Translate() override;

    bool importModel(int, const Handle(AIS_InteractiveContext) &);

    bool exportModel(int, const Handle(AIS_InteractiveContext) &);

    QString info() const;

protected:
    virtual Handle(TopTools_HSequenceOfShape) importModel(int, const QString &);

    virtual bool exportModel(int, const QString &,
                             const Handle(TopTools_HSequenceOfShape) &);

    virtual bool displayShSequence(const Handle(AIS_InteractiveContext) &,
                                   const Handle(TopTools_HSequenceOfShape) &);

    QString selectFileName(int, bool);

private:
    TranslateDlg *getDialog(int, bool);

    static Handle(TopTools_HSequenceOfShape) getShapes(const Handle(AIS_InteractiveContext) &);

    static Handle(TopTools_HSequenceOfShape) importBREP(const QString &);

    static Handle(TopTools_HSequenceOfShape) importIGES(const QString &);

    static Handle(TopTools_HSequenceOfShape) importSTEP(const QString &);

    static bool exportBREP(const QString &, const Handle(TopTools_HSequenceOfShape) &);

    static bool exportIGES(const QString &, const Handle(TopTools_HSequenceOfShape) &);

    bool exportSTEP(const QString &, const Handle(TopTools_HSequenceOfShape) &);

    bool exportSTL(const QString &, const Handle(TopTools_HSequenceOfShape) &);

    bool exportVRML(const QString &, const Handle(TopTools_HSequenceOfShape) &);

    static bool checkFacetedBrep(const Handle(TopTools_HSequenceOfShape) &);

protected:
    TranslateDlg *myDlg;
    QString myInfo;
};

#endif

