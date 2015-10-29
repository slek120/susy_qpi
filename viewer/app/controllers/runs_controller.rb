class RunsController < ApplicationController
  def index
    @runs = Run.order('date desc, time desc').page params[:page]
    @runs = @runs.where(t: params[:t]) unless params[:t].blank?
    @runs = @runs.where(mu: params[:mu]) unless params[:mu].blank?
    @runs = @runs.where(x0: params[:x0]) unless params[:x0].blank?
    @runs = @runs.where(epsf: params[:epsf]) unless params[:epsf].blank?
    @runs = @runs.where(V: params[:V]) unless params[:V].blank?
    @runs = @runs.where(Uc: params[:Uc]) unless params[:Uc].blank?
    @runs = @runs.where(Uf: params[:Uf]) unless params[:Uf].blank?
    @runs = @runs.where(omega: params[:omega]) unless params[:omega].blank?
    @runs = @runs.where(delta: params[:delta]) unless params[:delta].blank?

    case params[:order]
    when 'tasc'
      @runs = @runs.order('t asc')
    when 'tdesc'
      @runs = @runs.order('t desc')
    when 'muasc'
      @runs = @runs.order('mu asc')
    when 'mudesc'
      @runs = @runs.order('mu desc')
    when 'x0asc'
      @runs = @runs.order('x0 asc')
    when 'x0desc'
      @runs = @runs.order('x0 desc')
    when 'epsfasc'
      @runs = @runs.order('epsf asc')
    when 'epsfdesc'
      @runs = @runs.order('epsf desc')
    when 'Vasc'
      @runs = @runs.order('V asc')
    when 'Vdesc'
      @runs = @runs.order('V desc')
    when 'Ucasc'
      @runs = @runs.order('Uc asc')
    when 'Ucdesc'
      @runs = @runs.order('Uc desc')
    when 'Ufasc'
      @runs = @runs.order('Uf asc')
    when 'Ufdesc'
      @runs = @runs.order('Uf desc')
    when 'omegaasc'
      @runs = @runs.order('omega asc')
    when 'omegadesc'
      @runs = @runs.order('omega desc')
    when 'deltaasc'
      @runs = @runs.order('delta asc')
    when 'deltadesc'
      @runs = @runs.order('delta desc')
    end
  end
end
